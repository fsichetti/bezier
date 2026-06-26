#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVectorT<3, 3, 1>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 24);
	assert(out.size() == 4);
	std::array<Interval, 24> cp;
	for (int i = 0; i < 24; ++i) cp[i] = cpFP[i];
	I tmp_0 = -cp[4];
	I tmp_1 = cp[10] + tmp_0;
	I tmp_2 = -cp[2];
	I tmp_3 = cp[20] + tmp_2;
	I tmp_4 = cp[22] + tmp_0;
	I tmp_5 = -cp[2] + cp[8];
	I tmp_6 = cp[14] + tmp_2;
	I tmp_7 = cp[16] + tmp_0;
	I tmp_8 = (R(2, 3))*cp[0] + (R(1, 3))*cp[1];
	I tmp_9 = -R(2, 3)*cp[4] - R(1, 3)*cp[5];
	I tmp_10 = (R(2, 3))*cp[10] + (R(1, 3))*cp[11] + tmp_9;
	I tmp_11 = (R(2, 3))*cp[2];
	I tmp_12 = (R(1, 3))*cp[3];
	I tmp_13 = -tmp_11 - tmp_12;
	I tmp_14 = (R(2, 3))*cp[14] + (R(1, 3))*cp[15] + tmp_13;
	I tmp_15 = (R(2, 3))*cp[16] + (R(1, 3))*cp[17] + tmp_9;
	I tmp_16 = (R(2, 3))*cp[8] + (R(1, 3))*cp[9] - tmp_11 - tmp_12;
	I tmp_17 = (R(2, 3))*cp[20] + (R(1, 3))*cp[21] + tmp_13;
	I tmp_18 = (R(2, 3))*cp[22] + (R(1, 3))*cp[23] + tmp_9;
	I tmp_19 = (R(1, 3))*cp[0] + (R(2, 3))*cp[1];
	I tmp_20 = -R(1, 3)*cp[4] - R(2, 3)*cp[5];
	I tmp_21 = (R(1, 3))*cp[10] + (R(2, 3))*cp[11] + tmp_20;
	I tmp_22 = (R(2, 3))*cp[3];
	I tmp_23 = (R(1, 3))*cp[2];
	I tmp_24 = -tmp_22 - tmp_23;
	I tmp_25 = (R(1, 3))*cp[14] + (R(2, 3))*cp[15] + tmp_24;
	I tmp_26 = (R(1, 3))*cp[16] + (R(2, 3))*cp[17] + tmp_20;
	I tmp_27 = (R(1, 3))*cp[8] + (R(2, 3))*cp[9] - tmp_22 - tmp_23;
	I tmp_28 = (R(1, 3))*cp[20] + (R(2, 3))*cp[21] + tmp_24;
	I tmp_29 = (R(1, 3))*cp[22] + (R(2, 3))*cp[23] + tmp_20;
	I tmp_30 = -cp[1];
	I tmp_31 = -cp[5];
	I tmp_32 = cp[11] + tmp_31;
	I tmp_33 = -cp[3];
	I tmp_34 = cp[21] + tmp_33;
	I tmp_35 = cp[23] + tmp_31;
	I tmp_36 = -cp[3] + cp[9];
	I tmp_37 = cp[15] + tmp_33;
	I tmp_38 = cp[17] + tmp_31;
	out[0] = (-cp[0] + cp[12])*(tmp_1*tmp_3 - tmp_4*tmp_5) + (-cp[0] + cp[18])*(-tmp_1*tmp_6 + tmp_5*tmp_7) + (-cp[0] + cp[6])*(-tmp_3*tmp_7 + tmp_4*tmp_6);
	out[1] = (-tmp_10*tmp_14 + tmp_15*tmp_16)*((R(2, 3))*cp[18] + (R(1, 3))*cp[19] - tmp_8) + (tmp_10*tmp_17 - tmp_16*tmp_18)*((R(2, 3))*cp[12] + (R(1, 3))*cp[13] - tmp_8) + (tmp_14*tmp_18 - tmp_15*tmp_17)*((R(2, 3))*cp[6] + (R(1, 3))*cp[7] - tmp_8);
	out[2] = (-tmp_21*tmp_25 + tmp_26*tmp_27)*((R(1, 3))*cp[18] + (R(2, 3))*cp[19] - tmp_19) + (tmp_21*tmp_28 - tmp_27*tmp_29)*((R(1, 3))*cp[12] + (R(2, 3))*cp[13] - tmp_19) + (tmp_25*tmp_29 - tmp_26*tmp_28)*((R(1, 3))*cp[6] + (R(2, 3))*cp[7] - tmp_19);
	out[3] = (cp[13] + tmp_30)*(tmp_32*tmp_34 - tmp_35*tmp_36) + (cp[19] + tmp_30)*(-tmp_32*tmp_37 + tmp_36*tmp_38) + (-cp[1] + cp[7])*(-tmp_34*tmp_38 + tmp_35*tmp_37);
}
}
#undef R
#undef I