#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<3, 3, 1>(
	const span<const fp_t> xFP,
	const span<const Interval> lagVec
) {
	assert(xFP.size() == 3);
	assert(lagVec.size() == 1);
	std::array<Interval, 3> x;
	for (int i = 0; i < 3; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	acc += lagVec[0] * Interval(1);
	return acc;
}}
#undef R
#undef I