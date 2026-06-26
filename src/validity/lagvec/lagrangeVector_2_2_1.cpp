#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<2, 2, 1>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 6);
	assert(out.size() == 1);
	std::array<Interval, 6> cp;
	for (int i = 0; i < 6; ++i) cp[i] = cpFP[i];
	out[0] = (-cp[0] + cp[2])*(-cp[1] + cp[5]) - (-cp[0] + cp[4])*(-cp[1] + cp[3]);
}
}
#undef R
#undef I