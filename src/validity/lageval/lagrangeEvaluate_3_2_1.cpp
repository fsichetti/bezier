#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<3, 2, 1>(
	const span<const fp_t> xFP,
	const span<const Interval> lagVec
) {
	assert(xFP.size() == 3);
	assert(lagVec.size() == 9);
	std::array<Interval, 3> x;
	for (int i = 0; i < 3; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = powi(x[2], 2);
	I tmp_1 = -2*tmp_0;
	I tmp_2 = 3*x[2];
	I tmp_3 = 2*tmp_0;
	I tmp_4 = tmp_3*x[1];
	I tmp_5 = -tmp_2*x[1] + tmp_4 + x[1];
	I tmp_6 = tmp_3*x[0];
	I tmp_7 = -tmp_2*x[0] + tmp_6 + x[0];
	I tmp_8 = 4*tmp_0;
	I tmp_9 = 4*x[2];
	I tmp_10 = tmp_8*x[1] - tmp_9*x[1];
	I tmp_11 = tmp_8*x[0] - tmp_9*x[0];
	I tmp_12 = tmp_4 - x[1]*x[2];
	I tmp_13 = tmp_6 - x[0]*x[2];
	acc += lagVec[0] * (-tmp_1 - tmp_2 - tmp_5 - tmp_7 + 1);
	acc += lagVec[1] * (tmp_10 + tmp_11 - tmp_8 + tmp_9);
	acc += lagVec[2] * (-tmp_1 - tmp_12 - tmp_13 - x[2]);
	acc += lagVec[3] * (tmp_5);
	acc += lagVec[4] * (-tmp_10);
	acc += lagVec[5] * (tmp_12);
	acc += lagVec[6] * (tmp_7);
	acc += lagVec[7] * (-tmp_11);
	acc += lagVec[8] * (tmp_13);
	return acc;
}}
#undef R
#undef I