#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVectorT<1, 1, 1>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 4);
	assert(out.size() == 2);
	std::array<Interval, 4> cp;
	for (int i = 0; i < 4; ++i) cp[i] = cpFP[i];
	out[0] = -cp[0] + cp[2];
	out[1] = -cp[1] + cp[3];
}
}
#undef R
#undef I