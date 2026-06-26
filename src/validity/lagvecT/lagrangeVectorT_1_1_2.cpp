#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVectorT<1, 1, 2>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 6);
	assert(out.size() == 4);
	std::array<Interval, 6> cp;
	for (int i = 0; i < 6; ++i) cp[i] = cpFP[i];
	I tmp_0 = -4*cp[2];
	I tmp_1 = -4*cp[3];
	out[0] = -3*cp[0] - cp[4] - tmp_0;
	out[1] = -3*cp[1] - cp[5] - tmp_1;
	out[2] = cp[0] + 3*cp[4] + tmp_0;
	out[3] = cp[1] + 3*cp[5] + tmp_1;
}
}
#undef R
#undef I