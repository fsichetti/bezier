#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<1, 1, 4>(
	const span<const fp_t> xFP,
	const span<const Interval> lagVec
) {
	assert(xFP.size() == 1);
	assert(lagVec.size() == 4);
	std::array<Interval, 1> x;
	for (int i = 0; i < 1; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = powi(x[0], 2);
	I tmp_1 = powi(x[0], 3);
	I tmp_2 = (R(9, 2))*tmp_1;
	I tmp_3 = (R(27, 2))*tmp_1;
	acc += lagVec[0] * (9*tmp_0 - tmp_2 - R(11, 2)*x[0] + 1);
	acc += lagVec[1] * (-R(45, 2)*tmp_0 + tmp_3 + 9*x[0]);
	acc += lagVec[2] * (18*tmp_0 - tmp_3 - R(9, 2)*x[0]);
	acc += lagVec[3] * (-R(9, 2)*tmp_0 + tmp_2 + x[0]);
	return acc;
}}
#undef R
#undef I