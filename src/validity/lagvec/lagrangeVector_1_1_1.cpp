#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<1, 1, 1>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 2);
	assert(out.size() == 1);
	std::array<Interval, 2> cp;
	for (int i = 0; i < 2; ++i) cp[i] = cpFP[i];
	out[0] = -cp[0] + cp[1];
}
}
#undef R
#undef I