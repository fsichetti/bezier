#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<1, 1, 2>(
	const span<const fp_t> xFP,
	const span<const Interval> lagVec
) {
	assert(xFP.size() == 1);
	assert(lagVec.size() == 2);
	std::array<Interval, 1> x;
	for (int i = 0; i < 1; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	acc += lagVec[0] * (1 - x[0]);
	acc += lagVec[1] * (x[0]);
	return acc;
}}
#undef R
#undef I