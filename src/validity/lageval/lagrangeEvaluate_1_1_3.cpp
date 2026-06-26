#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<1, 1, 3>(
	const span<const fp_t> xFP,
	const span<const Interval> lagVec
) {
	assert(xFP.size() == 1);
	assert(lagVec.size() == 3);
	std::array<Interval, 1> x;
	for (int i = 0; i < 1; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = powi(x[0], 2);
	I tmp_1 = 2*tmp_0;
	acc += lagVec[0] * (tmp_1 - 3*x[0] + 1);
	acc += lagVec[1] * (-4*tmp_0 + 4*x[0]);
	acc += lagVec[2] * (tmp_1 - x[0]);
	return acc;
}}
#undef R
#undef I