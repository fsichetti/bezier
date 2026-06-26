#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<1, 1, 2>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 3);
	assert(out.size() == 2);
	std::array<Interval, 3> cp;
	for (int i = 0; i < 3; ++i) cp[i] = cpFP[i];
	I tmp_0 = -4*cp[1];
	out[0] = -3*cp[0] - cp[2] - tmp_0;
	out[1] = cp[0] + 3*cp[2] + tmp_0;
}
}
#undef R
#undef I