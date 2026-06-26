#pragma once
#include "Validator.hpp"
#include <utils/par_for.hpp>
#include <atomic>
namespace element_validity {
template<int n, int s, int p>
class StaticValidator : public Validator {
	private:
	static constexpr int subdomains = 1U << n;
	std::array<Matrix<Interval>, subdomains> matQ;

	Validity isValidElement(
		span<const fp_t> cp,
		std::vector<int> *adaptiveHierarchy = nullptr,
		Info *info = nullptr
	) const;
	Validity isValidMesh(
		span<const fp_t> cp,
		std::vector<int> *adaptiveHierarchy = nullptr,
		int *invalidElemID = nullptr,
		std::vector<int> *invalidList = nullptr,
		bool globalTest = true
	) const;

	Subdomain split(const Subdomain &src, int q) const;

	public:
	StaticValidator(int nThreads = 1) : Validator(nThreads) {
		initMatrices<n, s, p>(matL2B, matQ);
		cornerIndices<n, s, p>(interpIndices);
	}
	
	Validity isValid(
		span<const fp_t> cp,
		std::vector<int> *adaptiveHierarchy = nullptr,
		int *invalidElemID = nullptr,
		std::vector<int> *invalidList = nullptr,
		bool globalTest = true,
		Info *info = nullptr
	) const;

	Validity isValidAtTime(
		span<const fp_t> cp,
		fp_t time,
		std::vector<int> *adaptiveHierarchy = nullptr,
		int *invalidElemID = nullptr,
		std::vector<int> *invalidList = nullptr,
		bool globalTest = true,
		Info *info = nullptr
	) const;

	#ifdef EIGEN_INTERFACE
	Validity isValid(
		const Eigen::MatrixXd& cp,
		std::vector<int> *adaptiveHierarchy = nullptr,
		int *invalidElemID = nullptr,
		std::vector<int> *invalidList = nullptr,
		bool globalTest = true,
		Info *info = nullptr
	) const {
		return isValid(
			convertEigenMatrix(cp),
			adaptiveHierarchy,
			invalidElemID,
			invalidList,
			info
		);
	}
	#endif
};

//------------------------------------------------------------------------------

template<int n, int s, int p>
Validator::Subdomain StaticValidator<n,s,p>::split(
	const Validator::Subdomain &src, int q
) const {
	Subdomain res;
	matQ[q].mult(src.B, res.B);
	res.incl = minclusion(res.B);
	src.copySequence(res.qSequence);
	res.qSequence.push_back(q);
	return res;
}

//------------------------------------------------------------------------------

template<int n, int s, int p>
Validity StaticValidator<n,s,p>::isValid(
	span<const fp_t> cp,
	std::vector<int> *adaptiveHierarchy,
	int *invalidElemID,
	std::vector<int> *invalidList,
	bool globalTest,
	Info *info
) const {
	constexpr int numCoordsPerElem = nControlGeoMap(n,s,p) * n;
	const int numEl = cp.size() / numCoordsPerElem;
	if (numEl == 1) return isValidElement(cp, adaptiveHierarchy, info);
	else return isValidMesh(cp, adaptiveHierarchy, invalidElemID, invalidList, globalTest);
}

template<int n, int s, int p>
Validity StaticValidator<n,s,p>::isValidAtTime(
	span<const fp_t> cp,
	fp_t time,
	std::vector<int> *adaptiveHierarchy,
	int *invalidElemID,
	std::vector<int> *invalidList,
	bool globalTest,
	Info *info
) const {
	std::vector<fp_t> frame;
	assert(cp.size() % 2 == 0);
	frame.reserve(cp.size() / 2);
	const fp_t time1m = 1.-time;
	for (int i = 0; i < cp.size(); i+=2)
		frame.push_back(time1m*cp[i] + time*cp[i+1]);
	return isValid(frame, adaptiveHierarchy, invalidElemID, invalidList, globalTest, info);
}

//------------------------------------------------------------------------------

template<int n, int s, int p>
Validity StaticValidator<n, s, p>::isValidElement(
	span<const fp_t> cp,
	std::vector<int> *hierarchy,
	Info *info
) const {
	// Compute Lagrange coefficients
	std::vector<Interval> vL(nControlJacobian(n,s,p,false));
	lagrangeVector<n, s, p>(cp, vL);

	// Check for early termination
	for (const Interval &l : vL) {
		if (l < 0) {
			if (info) info->status = Info::Status::reachedTarget;
			return Validity::invalid;
		}
	}


	// Create initial subdomain
	Subdomain sd0;
	matL2B.mult(vL, sd0.B);
	sd0.incl = minclusion(sd0.B);
	
	// Initialize queue
	std::priority_queue<Subdomain> queue;
	queue.push(sd0);

	// Initialize auxiliary variables
	int reachedDepth = 0;
	const bool maxIterCheck = (maxSubdiv > 0);
	bool foundInvalid = false;
	bool gaveUp = false;

	if (hierarchy) hierarchy->clear();

	// Loop
	while(true) {
		if (queue.empty()) {
			if (info) info->status = Info::Status::completed;
			break;
		}
	
		// Get box from top of the queue
		const Subdomain dom = queue.top();
		queue.pop();

		reachedDepth = std::max(reachedDepth, dom.depth());

		// Check whether we need to give up
		if (maxIterCheck && (reachedDepth >= maxSubdiv)) {
			gaveUp = true;
			if (hierarchy && !foundInvalid) dom.copySequence(*hierarchy);
			if (info) info->status = Info::Status::maxDepth;
			break;
		}

		// Subdomain contains invalidity
		if (dom.incl <= epsilon) {
			foundInvalid = true;
			if (info) info->status = Info::Status::reachedTarget;
		}
		// Subdomain is undetermined and needs refinement
		else if (!(dom.incl > epsilon)) {
			// Split on all axes and push to queue
			for (int q=0; q<subdomains; ++q) queue.push(split(dom, q));
		}
	}

	if (info) info->spaceDepth = reachedDepth;
	return gaveUp ? Validity::uncertain : (foundInvalid ? Validity::invalid : Validity::valid);
}

//------------------------------------------------------------------------------

template<int n, int s, int p>
Validity StaticValidator<n, s, p>::isValidMesh(
	span<const fp_t> cp,
	std::vector<int> *adaptiveHierarchy,
	int *invalidElemID,
	std::vector<int> *invalidList,
	bool globalTest
) const {
	const int numCoordsPerElem = nControlGeoMap(n,s,p) * n;
	const int numEl = cp.size() / (numCoordsPerElem);
	std::vector<Validity> results(numEl);
	std::vector<fp_t> timings(numEl);
	std::vector<std::vector<int>> hierarchies(numEl);
	bool gaveUp = false;
	std::atomic<bool> foundInvalid = false;

	par_for(numEl, nThreads, [&](int start, int end, int thread_id) {
		for (int e = start; e < end; ++e) {
			if (foundInvalid) {
				timings.at(e) = 0;
				continue;
			}
			span<const fp_t> element(
				cp.data() + numCoordsPerElem * e, numCoordsPerElem);
			Timer timer;
			timer.start();
			Validity v = isValidElement(element, &hierarchies.at(e));
			timer.stop();
			results.at(e) = v;
			if (v == Validity::invalid) foundInvalid = true;

			timings.at(e) = timer.read<std::chrono::microseconds>();
		}
	});

	for (int e = 0; e < numEl; e++)
		if (results.at(e) == Validity::uncertain)
			gaveUp = true;

	if (invalidList) {
		invalidList->clear();
		for (int e = 0; e < numEl; e++) {
			if (results.at(e) == Validity::invalid || results.at(e) == Validity::uncertain) {
				invalidList->push_back(e);
			}
		}
	}

	if (foundInvalid) {
		int invalidIndex = 0;
		while (results.at(invalidIndex) != Validity::invalid) ++invalidIndex;
		if (invalidElemID) *invalidElemID = invalidIndex;
		if (adaptiveHierarchy) *adaptiveHierarchy = std::move(hierarchies.at(invalidIndex));
		return Validity::invalid;
	}
	else if (gaveUp) {
		std::vector<int> depthOfSequence(numEl);
		for (int j = 0; j < numEl; ++j)
			depthOfSequence.at(j) = (results.at(j) == Validity::valid) ?
				-1 : hierarchies.at(j).size();
		const int i = std::distance(depthOfSequence.cbegin(),
			std::max_element(
				depthOfSequence.cbegin(),
				depthOfSequence.cend()
			)
		);
		if (invalidElemID) *invalidElemID = i;
		if (adaptiveHierarchy) *adaptiveHierarchy = std::move(hierarchies.at(i));
		return Validity::uncertain;
	}
	return Validity::valid;
}

}
