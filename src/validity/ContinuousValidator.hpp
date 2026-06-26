#pragma once
#include "Validator.hpp"
#include <utils/par_for.hpp>
#include <atomic>
namespace element_validity {

template<int n, int s, int p>
class ContinuousValidator : public Validator {
	private:
	static constexpr int subdomains = 2U << n;
	std::array<Matrix<Interval>, 2> matT;
	std::array<Matrix<Interval>, subdomains> matQ;

	fp_t precision = .1;

	fp_t maxTimeStepElement(
		span<const fp_t> cp,
		std::vector<int> *adaptiveHierarchy = nullptr,
		fp_t *timeOfInversion = nullptr,
		fp_t earlyStop = 1.,
		Info *info = nullptr
	) const;

	fp_t maxTimeStepMesh(
		span<const fp_t> cp,
		std::vector<int> *adaptiveHierarchy = nullptr,
		int *invalidElemID = nullptr,
		fp_t *timeOfInversion = nullptr
	) const;
	
	Subdomain split(const Subdomain &src, int q) const;

	public:
	ContinuousValidator(int nThreads = 1) : Validator(nThreads) {
		initMatricesT<n, s, p>(matL2B, matT, matQ);
		cornerIndicesT<n, s, p>(interpIndices);
	}
	fp_t maxTimeStep(
		span<const fp_t> cp,
		std::vector<int> *adaptiveHierarchy = nullptr,
		int *invalidElemID = nullptr,
		fp_t *timeOfInversion = nullptr,
		Info *info = nullptr
	) const;

	#ifdef EIGEN_INTERFACE
	fp_t maxTimeStep(
		const Eigen::MatrixXd& cp0,
		const Eigen::MatrixXd& cp1,
		std::vector<int> *adaptiveHierarchy = nullptr,
		int *invalidElemID = nullptr,
		fp_t *timeOfInversion = nullptr,
		Info *info = nullptr
	) const {
		return maxTimeStep(
			convertEigenMatrix(cp0, cp1),
			adaptiveHierarchy,
			invalidElemID,
			timeOfInversion,
			info
		);
	}
	#endif

	void setPrecisionTarget(fp_t t) { precision = t; }

	private:
	Subdomain splitTime(const Subdomain &src, bool t) const;
};

//------------------------------------------------------------------------------

template<int n, int s, int p>
Validator::Subdomain ContinuousValidator<n, s, p>::split(
	const Subdomain &src, int q
) const {
	Subdomain res;
	matQ[q].mult(src.B, res.B);
	res.time = (q >> n) ?
		RealInterval(src.time.middle(), src.time.upper()) :
		RealInterval(src.time.lower(), src.time.middle());
	res.incl = minclusion(res.B);
	src.copySequence(res.qSequence);
	res.qSequence.push_back(q % (1U << n));
	return res;
}

template<int n, int s, int p>
Validator::Subdomain ContinuousValidator<n, s, p>::splitTime(
	const Subdomain &src, bool t
) const {
	Subdomain res;
	(t ? matT[1] : matT[0]).mult(src.B, res.B);
	res.time = t ?
		RealInterval(src.time.middle(), src.time.upper()) :
		RealInterval(src.time.lower(), src.time.middle());
	res.incl = minclusion(res.B);
	src.copySequence(res.qSequence);
	return res;
}

//------------------------------------------------------------------------------

template<int n, int s, int p>
fp_t ContinuousValidator<n,s,p>::maxTimeStep(
	span<const fp_t> cp,
	std::vector<int> *adaptiveHierarchy,
	int *invalidElemID,
	fp_t *timeOfInversion,
	Info *info
) const {
	const int numCoordsPerElem = nControlGeoMap(n,s,p) * 2 * n;
	const int numEl = cp.size() / (numCoordsPerElem);
	if (numEl == 1) return maxTimeStepElement(
		cp, adaptiveHierarchy, timeOfInversion, 1., info);
	else return maxTimeStepMesh(
		cp, adaptiveHierarchy, invalidElemID, timeOfInversion);
}

//------------------------------------------------------------------------------

template<int n, int s, int p>
fp_t ContinuousValidator<n, s, p>::maxTimeStepElement(
	span<const fp_t> cp,
	std::vector<int> *hierarchy,
	fp_t *timeOfInversion,
	fp_t earlyStop,
	Info *info
) const {
	assert(precision > 0 && precision <= 1);

	// Compute Lagrange coefficients
	std::vector<Interval> vL(nControlJacobian(n,s,p,true));
	lagrangeVectorT<n, s, p>(cp, vL);

	// Create initial subdomain
	Subdomain sd0;
	matL2B.mult(vL, sd0.B);
	sd0.time = RealInterval(0,1);
	sd0.incl = minclusion(sd0.B);
	
	// Initialize queue
	std::priority_queue<Subdomain> queue;
	queue.push(sd0);

	// Initialize auxiliary variables
	int reachedDepthS = 0;
	int reachedDepthT = 0;
	const bool maxIterCheck = (maxSubdiv > 0);
	bool foundInvalid = false;
	fp_t tmin = 0.;
	fp_t tmax = 1.;

	if (hierarchy) hierarchy->clear();

	// Loop
	while(true) {
		if (queue.empty()) {
			tmin = 1;
			if (info) info->status = Info::Status::completed;
			break;
		}
		// Check whether we reached precision
		if (foundInvalid && tmin > 0. && tmax - tmin < precision) {
			if (info) info->status = Info::Status::reachedTarget;
			break;
		}

		// Check whether we satisfy early termination
		if (tmin >= earlyStop) {
			if (info) info->status = Info::Status::globalCondition;
			break;
		}

		// Get box from top of the queue
		const Subdomain dom = queue.top();
		queue.pop();
		// for (const auto u : dom.qSequence) std::cerr << u << "-";
		// std::cerr << dom.time << "\t" << dom.incl << std::endl;

		if (dom.time.width() == 0) {
			if (info) info->status = Info::Status::noSplit;
			break;
		}
		assert(dom.time.lower() < tmax);
		assert(dom.time.lower() >= tmin);
		reachedDepthS = std::max(reachedDepthS, dom.depth());
		reachedDepthT = std::max(reachedDepthT, dom.timeDepth());

		// Check whether we need to give up
		if (maxIterCheck && (reachedDepthS >= maxSubdiv)) {
			if (hierarchy && !foundInvalid) dom.copySequence(*hierarchy);
			if (info) info->status = Info::Status::maxDepth;
			break;
		}

		// Update tmin
		tmin = dom.time.lower();
		if (!(dom.incl > epsilon)) {
			if constexpr (p == 1 && n == s) {
				// Subdomain contains invalidity
				if (dom.incl <= epsilon) {
					foundInvalid = true;
					if (tmax > dom.time.upper()) {
						tmax = dom.time.upper();				// update tmax
						if (hierarchy) dom.copySequence(*hierarchy);
					}
				}
				// Subdomain contains invalidity or is uncertain
				queue.push(splitTime(dom, false));
				queue.push(splitTime(dom, true));
			}
			else {
				#ifdef EXPERIM_NAIVE_SUBDIV
					// Subdomain contains invalidity
					if (dom.incl <= epsilon) {
						foundInvalid = true;
						if (tmax > dom.time.upper()) {
							tmax = dom.time.upper();				// update tmax
							if (hierarchy) dom.copySequence(*hierarchy);
						}
					}
					// Split on all axes and push to queue
					for (int q=0; q<subdomains; ++q) queue.push(split(dom, q));
				#else
					// Subdomain contains invalidity
					if (dom.incl <= epsilon) {
						foundInvalid = true;
						if (tmax > dom.time.upper()) {
							tmax = dom.time.upper();				// update tmax
							if (hierarchy) dom.copySequence(*hierarchy);
						}
						// Split on time only and push to queue
						queue.push(splitTime(dom, false));
						queue.push(splitTime(dom, true));
					}
					// Subdomain is undetermined and needs refinement
					else {
						// Split on all axes and push to queue
						for (int q=0; q<subdomains; ++q)
							queue.push(split(dom, q));
					}
				#endif
			}
		}
	}

	if (info) info->spaceDepth = reachedDepthS;
	if (info) info->timeDepth = reachedDepthT;
	if (timeOfInversion) *timeOfInversion = foundInvalid ? tmax : 2.;
	return tmin;
}

//------------------------------------------------------------------------------

template<int n, int s, int p>
fp_t ContinuousValidator<n, s, p>::maxTimeStepMesh(
	span<const fp_t> cp,
	std::vector<int> *adaptiveHierarchy,
	int *invalidElemID,
	fp_t *toi
) const {
	const int numCoordsPerElem = nControlGeoMap(n,s,p) * 2 * n;
	const int numEl = cp.size() / (numCoordsPerElem);
	std::atomic<fp_t> minT = 1;
	std::vector<fp_t> timeOfInversion(numEl);
	std::vector<fp_t> timings(numEl);
	std::vector<std::vector<int>> hierarchies(numEl);

	par_for(numEl, nThreads, [&](int start, int end, int thread_id) {
		for (int e = start; e < end; ++e) {
			span<const fp_t> element(
				cp.data() + numCoordsPerElem * e, numCoordsPerElem);
			Timer timer;
			timer.start();
			fp_t t = maxTimeStepElement(
				element,
				&hierarchies.at(e),
				&timeOfInversion.at(e),
				minT
			);
			timer.stop();
			timings.at(e) = timer.read<std::chrono::microseconds>();
			if (t < minT)
				minT = t;
		}
	});

	bool foundInvalid = false;
	for (int e = 0; e < numEl; e++) {
		if (timeOfInversion.at(e) <= 1.) {
			foundInvalid = true;
			break;
		}
	}

	if (foundInvalid) {
		const auto m =
			std::min_element(timeOfInversion.cbegin(), timeOfInversion.cend()); 
		const int i = std::distance(timeOfInversion.cbegin(), m);
		if (toi) *toi = *m;
		if (invalidElemID) *invalidElemID = i;
		if (adaptiveHierarchy) *adaptiveHierarchy = std::move(hierarchies.at(i));
	}
	else {
		std::vector<fp_t> depthOfSequence(numEl);
		for (int j = 0; j < numEl; ++j)
			depthOfSequence.at(j) = hierarchies.at(j).size();
		const auto m =
			std::max_element(depthOfSequence.cbegin(), depthOfSequence.cend());
		const int i = std::distance(depthOfSequence.cbegin(), m);
		if (invalidElemID) *invalidElemID = i;
		if (adaptiveHierarchy) *adaptiveHierarchy = std::move(hierarchies.at(i));
	}
	return minT;
}


}