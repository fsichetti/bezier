#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVectorT<2, 2, 1>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 12);
	assert(out.size() == 3);
	std::array<Interval, 12> cp;
	for (int i = 0; i < 12; ++i) cp[i] = cpFP[i];
	I tmp_0 = (R(1, 2))*cp[0] + (R(1, 2))*cp[1];
	I tmp_1 = (R(1, 2))*cp[2];
	I tmp_2 = (R(1, 2))*cp[3];
	out[0] = (-cp[0] + cp[4])*(cp[10] - cp[2]) - (-cp[0] + cp[8])*(-cp[2] + cp[6]);
	out[1] = ((R(1, 2))*cp[4] + (R(1, 2))*cp[5] - tmp_0)*((R(1, 2))*cp[10] + (R(1, 2))*cp[11] - tmp_1 - tmp_2) - ((R(1, 2))*cp[8] + (R(1, 2))*cp[9] - tmp_0)*((R(1, 2))*cp[6] + (R(1, 2))*cp[7] - tmp_1 - tmp_2);
	out[2] = (cp[11] - cp[3])*(-cp[1] + cp[5]) - (-cp[1] + cp[9])*(-cp[3] + cp[7]);
}
}
#undef R
#undef I