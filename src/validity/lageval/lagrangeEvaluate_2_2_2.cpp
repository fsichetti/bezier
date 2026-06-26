#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<2, 2, 2>(
	const span<const fp_t> xFP,
	const span<const Interval> lagVec
) {
	assert(xFP.size() == 2);
	assert(lagVec.size() == 6);
	std::array<Interval, 2> x;
	for (int i = 0; i < 2; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = 4*x[0]*x[1];
	I tmp_1 = powi(x[0], 2);
	I tmp_2 = 2*tmp_1;
	I tmp_3 = powi(x[1], 2);
	I tmp_4 = 2*tmp_3;
	acc += lagVec[0] * (tmp_0 + tmp_2 + tmp_4 - 3*x[0] - 3*x[1] + Interval(1.));
	acc += lagVec[1] * (-tmp_0 - 4*tmp_3 + 4*x[1]);
	acc += lagVec[2] * (tmp_4 - x[1]);
	acc += lagVec[3] * (-tmp_0 - 4*tmp_1 + 4*x[0]);
	acc += lagVec[4] * (tmp_0);
	acc += lagVec[5] * (tmp_2 - x[0]);
	return acc;
}}
#undef R
#undef I