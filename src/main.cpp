#include "element_validity.hpp"
#include "utils/Timer.hpp"
#include "utils/Settings.hpp"

#include <memory>
#include <fstream>

#if(defined GMP_INTERFACE && defined HDF5_INTERFACE)
#include <H5Cpp.h>

constexpr char SEP = ',';

template<
    int n, 
    int s, 
    int p
>
element_validity::fp_t processData(
    element_validity::Settings args,
    int nNodesPerElem,
    int numberOfElements,
    std::vector<element_validity::fp_t> nodes,
    std::ostream *out = nullptr
) {
    using namespace element_validity;
    Timer timer;

    const int nCoordPerElem = nNodesPerElem*n*2;
    const int nElements =
        args.numElem == 0 ? numberOfElements : args.numElem;
    const int lastElem = args.firstElem + nElements;
    fp_t minT = 1;

    ContinuousValidator<n, s, p> cChecker(args.numThreads);
    cChecker.setPrecisionTarget(args.precision);
    cChecker.setMaxSubdiv(args.maxIterations);
    StaticValidator<n, s, p> preChecker(args.numThreads);
    preChecker.setMaxSubdiv(args.preCheckMaxIter);
    StaticValidator<n, s, p> sChecker(args.numThreads);
    sChecker.setMaxSubdiv(args.maxIterations);
    
    if (!args.globalQuery) {
        if (out)
            *out << "ID" << SEP
                << "max_time_step" << SEP
                << "time_of_inversion" << SEP
                << "space_depth" << SEP
                << "time_depth" << SEP
                << "microseconds" << SEP
                << "hierarchy" << SEP
                << "description" << std::endl;

        // Continuous check
        for (int e=args.firstElem; e<lastElem; ++e) {
            const int elemOffset = e*nCoordPerElem;
            span<fp_t> element(nodes.data() + elemOffset, nCoordPerElem);
            std::vector<int> h;
            Validator::Info info;
            if (args.preCheck) {
                const Validity v0 = preChecker.isValidAtTime(element, 0);
                if (v0 != Validity::valid) {
                    if (v0 == Validity::uncertain)
                        std::cerr <<
                            "Warning: static check could not determine whether \
                            element " << e << " is valid at t=0" << std::endl;
                    continue;
                }
            }
            fp_t tInv = 2., mts = -1.;
            if (args.staticCheck) {
                timer.start();
                const Validity val = sChecker.isValidAtTime(
                    element, args.staticCheckTime, nullptr, nullptr, nullptr, true, &info);
                timer.stop();
                switch (val) {
                    case Validity::valid: mts = 1.; break;
                    case Validity::invalid: tInv = args.staticCheckTime; // don't break
                    case Validity::uncertain: mts = 0.; break;
                    default: throw std::runtime_error("Unknown validity value");
                }
            }
            else {
                timer.start();
                mts = cChecker.maxTimeStep(element, &h, nullptr, &tInv, &info);
                timer.stop();
            }
            minT = std::min(minT, mts);
            const double microseconds =
                static_cast<double>(timer.read<std::chrono::nanoseconds>()) / 1000;
            const std::string tInvStr = (tInv==2) ? "" : std::to_string(tInv);
            if (out) {
                *out << e << SEP;
                *out << fp_fmt << mts << SEP;
                *out << fp_fmt << tInvStr << SEP;
                *out << info.spaceDepth << SEP;
                *out << info.timeDepth << SEP;
                *out << microseconds << SEP;
                for (int u : h) *out << u << ' ';
                *out << SEP;
                *out << info.description() << std::endl;
            }
            timer.reset();
        }
    }
    else {
        std::vector<int> globalH;
        int invalidElemID;
        fp_t tInv;
        span<fp_t> elements(
            nodes.data() + args.firstElem*nCoordPerElem,
            nCoordPerElem*nElements
        );
        timer.start();
        minT = cChecker.maxTimeStep(elements, &globalH, &invalidElemID, &tInv);
        timer.stop();
        const double seconds =
            static_cast<double>(timer.read<std::chrono::milliseconds>()) / 1000;
        if (out) {
            *out << "Global info " << std::endl;
            *out << "Time in seconds: " << seconds << std::endl;
            *out << "Invalid element ID: " << invalidElemID << std::endl;
            *out << "Space subdivision hierarchy: ";
            for (int u : globalH) *out << u << ' ';
            *out << "Certified time of inversion: " << tInv << std::endl;
            *out << std::endl;
        }
    }
    return minT;
}

int main(int argc, char** argv) {
    using namespace element_validity;
    Interval::init();
    Settings args(argc, argv);
    if (args.abort) return 0;

    int dimension;
    int nNodesPerElem;
    int nElements;
    std::vector<fp_t> nodes;

    {
        // Read rational data from HDF5 dataset
        std::cout << "Reading file..." << std::flush;
        H5::H5File file(args.filePath, H5F_ACC_RDONLY);

        const auto H5int = H5::PredType::NATIVE_INT;
        file.openDataSet("Dimension").read(&dimension, H5int);
        file.openDataSet("NumberOfHighOrderNodes").read(&nNodesPerElem, H5int);
        file.openDataSet("NumberOfSimplices").read(&nElements, H5int);

        H5::DataSet dataset = file.openDataSet("Nodes");
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t nDataEntries = dataspace.getSimpleExtentNpoints();
        assert(4 * nElements * nNodesPerElem * dimension == nDataEntries);
        H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE);
        std::vector<char*> buffer(nDataEntries);
        dataset.read(buffer.data(), datatype, dataspace);
        nodes.reserve(nElements * nNodesPerElem * dimension * 2);
        for (hsize_t i=0; 2*i<nDataEntries; ++i) {
            Rational rat;
            rat.num() = buffer[2*i];
            rat.den() = buffer[2*i+1];
            nodes.push_back(static_cast<fp_t>(rat));
        }
        std::cout << " Done." << std::endl;
    }

    std::cout << "Processing..." << std::flush;
    std::unique_ptr<std::ostream> out;
    if (args.resultsPath.size() > 0)
        out = std::make_unique<std::ofstream>(args.resultsPath);
    std::ostream *const outptr = out ? out.get() : &std::cout;
    fp_t mts = 0;
    #define IFPROC(n, s, p) \
        if (dimension == n && nNodesPerElem == nControlGeoMap(n,s,p)) \
        mts = processData<n, s, p>(\
            args, nNodesPerElem, nElements, nodes, outptr);
    IFPROC(1, 1, 1)
    else IFPROC(1, 1, 2)
    else IFPROC(1, 1, 3)
    else IFPROC(1, 1, 4)
    else IFPROC(1, 1, 5)
    else IFPROC(2, 1, 1)
    else IFPROC(2, 1, 2)
    else IFPROC(2, 2, 1)
    else IFPROC(2, 2, 2)
    else IFPROC(2, 2, 3)
    else IFPROC(2, 2, 4)
    else IFPROC(3, 1, 1)
    else IFPROC(3, 2, 1)
    else IFPROC(3, 3, 1)
    else IFPROC(3, 3, 2)
    else IFPROC(3, 3, 3)
    #ifdef INSTANTIATE_3_3_4
    else IFPROC(3, 3, 4)
    #endif
    #undef IFPROC
    else throw std::invalid_argument("Not implemented");
    std::cout << " Done." << std::endl;
    std::cout << "Maximum valid time step: " << mts << std::endl;

    return 0;
}
#else

#ifdef _MSC_VER
    #pragma message("HDF5 or GMP interface disabled, ignoring main")
#elif defined(__GNUC__) || defined(__clang__)
    #warning "HDF5 or GMP interface disabled, ignoring main"
#endif

int main(int argc, char** argv) { return 0; }
#endif
