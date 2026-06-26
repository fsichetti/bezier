#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<2, 1, 1>(
	const span<const fp_t> xFP,
	const span<const Interval> lagVec
) {
	assert(xFP.size() == 2);
	assert(lagVec.size() == 4);
	std::array<Interval, 2> x;
	for (int i = 0; i < 2; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = -x[0];
	I tmp_1 = x[0]*x[1];
	I tmp_2 = tmp_1 - x[1];
	acc += lagVec[0] * (tmp_0 + tmp_2 + 1);
	acc += lagVec[1] * (-tmp_2);
	acc += lagVec[2] * (-tmp_0 - tmp_1);
	acc += lagVec[3] * (tmp_1);
	return acc;
}}
#undef R
#undef I