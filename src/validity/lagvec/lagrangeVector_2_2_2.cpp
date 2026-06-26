#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<2, 2, 2>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 12);
	assert(out.size() == 6);
	std::array<Interval, 12> cp;
	for (int i = 0; i < 12; ++i) cp[i] = cpFP[i];
	I tmp_0 = 3*cp[0];
	I tmp_1 = -4*cp[10];
	I tmp_2 = 3*cp[1];
	I tmp_3 = -4*cp[7];
	I tmp_4 = -4*cp[6];
	I tmp_5 = 4*cp[11];
	I tmp_6 = -tmp_5;
	I tmp_7 = cp[0] - cp[4];
	I tmp_8 = 2*cp[11];
	I tmp_9 = cp[1] - 2*cp[9];
	I tmp_10 = 2*cp[7];
	I tmp_11 = -tmp_10;
	I tmp_12 = cp[3] + tmp_11;
	I tmp_13 = cp[1] - cp[5];
	I tmp_14 = 2*cp[10];
	I tmp_15 = cp[0] - 2*cp[8];
	I tmp_16 = 2*cp[6];
	I tmp_17 = -tmp_16;
	I tmp_18 = cp[2] + tmp_17;
	I tmp_19 = 4*cp[8];
	I tmp_20 = cp[0] - cp[2];
	I tmp_21 = -tmp_8;
	I tmp_22 = cp[5] + tmp_21;
	I tmp_23 = -tmp_14;
	I tmp_24 = cp[4] + tmp_23;
	I tmp_25 = cp[0] + 2*cp[8];
	I tmp_26 = cp[1] + 2*cp[9];
	out[0] = (-cp[2] - tmp_0 - tmp_4)*(-cp[5] - tmp_2 - tmp_6) - (-cp[3] - tmp_2 - tmp_3)*(-cp[4] - tmp_0 - tmp_1);
	out[1] = -tmp_13*(-tmp_14 - tmp_15 - tmp_18) + tmp_7*(-tmp_12 - tmp_8 - tmp_9);
	out[2] = -(cp[0] + 3*cp[4] + tmp_1)*(cp[1] - cp[3] + 4*cp[9] - tmp_5) + (cp[1] + 3*cp[5] + tmp_6)*(tmp_1 + tmp_19 + tmp_20);
	out[3] = -tmp_20*(-tmp_10 - tmp_22 - tmp_9) - (-cp[1] + cp[3])*(-tmp_15 - tmp_16 - tmp_24);
	out[4] = (tmp_11 + tmp_22 + tmp_26)*(tmp_18 + tmp_23 + tmp_25) - (tmp_12 + tmp_21 + tmp_26)*(tmp_17 + tmp_24 + tmp_25);
	out[5] = (cp[0] + 3*cp[2] + tmp_4)*(4*cp[9] + tmp_13 + tmp_3) - (cp[1] + 3*cp[3] + tmp_3)*(tmp_19 + tmp_4 + tmp_7);
}
}
#undef R
#undef I