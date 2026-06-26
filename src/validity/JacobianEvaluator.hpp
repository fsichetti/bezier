#pragma once
#include "numbers/Interval.hpp"
#include "utils/combinatorics.hpp"
#include "lagrangeVector.hpp"
#include "lagrangeEvaluate.hpp"

namespace element_validity {
template<int n, int s, int p>
class JacobianEvaluator {
	private:
	std::vector<Interval> coeffs;
	static constexpr int numCoordsPerElem = nControlGeoMap(n,s,p) * n;
	static constexpr int numLagCoefPerElem = nControlJacobian(n,s,p,false);
	const int numEl = 0;

	public:
	JacobianEvaluator(const span<const fp_t>);
	fp_t eval(span<const fp_t>, int element=0) const;
	fp_t eval(std::initializer_list<fp_t> il, int element=0) const {
		std::vector<fp_t> v(il);
		return eval(v, element);
	}

	#ifdef EIGEN_INTERFACE
	JacobianEvaluator(const Eigen::MatrixXd& cp) :
		JacobianEvaluator(convertEigenMatrix(cp)) {}
	Eigen::VectorXd eval(const Eigen::MatrixXd& uv, int element=0) const {
		std::vector<double> uvs = convertEigenMatrix(uv);
		const int numPts = uvs.size() / n;
		Eigen::VectorXd res(numPts);
		for (int i=0; i<numPts; ++i) {
			span<fp_t> pt(uvs.data() + i*n, n);
			res[i] = eval(pt);
		}
		return res;
	}
	#endif
};

template<int n, int s, int p>
JacobianEvaluator<n, s, p>::JacobianEvaluator(const span<const fp_t> cp) :
	numEl(cp.size() / (numCoordsPerElem)
) {
	assert(cp.size() == numEl * numCoordsPerElem);
	coeffs.resize(numEl * numLagCoefPerElem);
	for (int i=0; i < numEl; ++i) {
		span<const fp_t> in =
			cp.subspan(i * numCoordsPerElem, numCoordsPerElem);
		// Save Lagrange coefficients to coeffs vector
		span<Interval> out(
			coeffs.data() + i * numLagCoefPerElem, numLagCoefPerElem);
		lagrangeVector<n, s, p>(in, out);
	}
}

template<int n, int s, int p>
fp_t JacobianEvaluator<n, s, p>::eval(span<const fp_t> x, int i) const {
	assert(x.size() == n);
	span<const Interval> lag(
		coeffs.data() + i * numLagCoefPerElem, numLagCoefPerElem);
	Interval res = lagrangeEvaluate<n,s,p>(x, lag);
	return res.middle();
}
}