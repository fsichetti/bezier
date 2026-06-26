#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<1, 1, 5>(
	const span<const fp_t> xFP,
	const span<const Interval> lagVec
) {
	assert(xFP.size() == 1);
	assert(lagVec.size() == 5);
	std::array<Interval, 1> x;
	for (int i = 0; i < 1; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = powi(x[0], 2);
	I tmp_1 = powi(x[0], 3);
	I tmp_2 = powi(x[0], 4);
	I tmp_3 = (R(32, 3))*tmp_2;
	I tmp_4 = (R(128, 3))*tmp_2;
	acc += lagVec[0] * ((R(70, 3))*tmp_0 - R(80, 3)*tmp_1 + tmp_3 - R(25, 3)*x[0] + 1);
	acc += lagVec[1] * (-R(208, 3)*tmp_0 + 96*tmp_1 - tmp_4 + 16*x[0]);
	acc += lagVec[2] * (76*tmp_0 - 128*tmp_1 + 64*tmp_2 - 12*x[0]);
	acc += lagVec[3] * (-R(112, 3)*tmp_0 + (R(224, 3))*tmp_1 - tmp_4 + (R(16, 3))*x[0]);
	acc += lagVec[4] * ((R(22, 3))*tmp_0 - 16*tmp_1 + tmp_3 - x[0]);
	return acc;
}}
#undef R
#undef I