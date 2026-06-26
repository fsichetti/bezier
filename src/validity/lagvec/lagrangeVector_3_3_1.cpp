#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<3, 3, 1>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 12);
	assert(out.size() == 1);
	std::array<Interval, 12> cp;
	for (int i = 0; i < 12; ++i) cp[i] = cpFP[i];
	I tmp_0 = cp[10] - cp[1];
	I tmp_1 = -cp[2] + cp[8];
	I tmp_2 = cp[11] - cp[2];
	I tmp_3 = -cp[1] + cp[7];
	I tmp_4 = -cp[2] + cp[5];
	I tmp_5 = -cp[1] + cp[4];
	out[0] = (-cp[0] + cp[3])*(-tmp_0*tmp_1 + tmp_2*tmp_3) + (-cp[0] + cp[6])*(tmp_0*tmp_4 - tmp_2*tmp_5) + (-cp[0] + cp[9])*(tmp_1*tmp_5 - tmp_3*tmp_4);
}
}
#undef R
#undef I