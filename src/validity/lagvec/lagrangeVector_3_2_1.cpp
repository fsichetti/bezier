#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<3, 2, 1>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 18);
	assert(out.size() == 9);
	std::array<Interval, 18> cp;
	for (int i = 0; i < 18; ++i) cp[i] = cpFP[i];
	I tmp_0 = -cp[3];
	I tmp_1 = -cp[0] - tmp_0;
	I tmp_2 = cp[10] - cp[1];
	I tmp_3 = -cp[8];
	I tmp_4 = -cp[2] - tmp_3;
	I tmp_5 = cp[11] - cp[2];
	I tmp_6 = -cp[7];
	I tmp_7 = -cp[1] - tmp_6;
	I tmp_8 = -cp[6];
	I tmp_9 = -cp[0] - tmp_8;
	I tmp_10 = -cp[5];
	I tmp_11 = -cp[2] - tmp_10;
	I tmp_12 = -cp[4];
	I tmp_13 = -cp[1] - tmp_12;
	I tmp_14 = -cp[9];
	I tmp_15 = -cp[0] - tmp_14;
	I tmp_16 = -tmp_11*tmp_7 + tmp_13*tmp_4;
	I tmp_17 = (R(1, 2))*cp[0] + (R(1, 2))*cp[9];
	I tmp_18 = (R(1, 2))*cp[15] + (R(1, 2))*cp[6] - tmp_17;
	I tmp_19 = (R(1, 2))*cp[11] + (R(1, 2))*cp[2];
	I tmp_20 = (R(1, 2))*cp[14] + (R(1, 2))*cp[5] - tmp_19;
	I tmp_21 = (R(1, 2))*cp[10] + (R(1, 2))*cp[1];
	I tmp_22 = (R(1, 2))*cp[13] + (R(1, 2))*cp[4] - tmp_21;
	I tmp_23 = (R(1, 2))*cp[12] + (R(1, 2))*cp[3] - tmp_17;
	I tmp_24 = (R(1, 2))*cp[17] + (R(1, 2))*cp[8] - tmp_19;
	I tmp_25 = (R(1, 2))*cp[16] + (R(1, 2))*cp[7] - tmp_21;
	I tmp_26 = -tmp_20*tmp_25 + tmp_22*tmp_24;
	I tmp_27 = -cp[10] + cp[13];
	I tmp_28 = -cp[11] + cp[17];
	I tmp_29 = -cp[10] + cp[16];
	I tmp_30 = -cp[11] + cp[14];
	I tmp_31 = tmp_27*tmp_28 - tmp_29*tmp_30;
	I tmp_32 = cp[12] + tmp_14;
	I tmp_33 = cp[15] + tmp_14;
	I tmp_34 = cp[16] + tmp_6;
	I tmp_35 = cp[17] + tmp_3;
	I tmp_36 = cp[15] + tmp_8;
	I tmp_37 = cp[13] + tmp_12;
	I tmp_38 = cp[14] + tmp_10;
	I tmp_39 = cp[12] + tmp_0;
	out[0] = tmp_1*(-tmp_2*tmp_4 + tmp_5*tmp_7) + tmp_15*tmp_16 + tmp_9*(tmp_11*tmp_2 - tmp_13*tmp_5);
	out[1] = tmp_15*tmp_26 + tmp_18*(tmp_2*tmp_20 - tmp_22*tmp_5) + tmp_23*(-tmp_2*tmp_24 + tmp_25*tmp_5);
	out[2] = tmp_15*tmp_31 + tmp_32*(-tmp_2*tmp_28 + tmp_29*tmp_5) + tmp_33*(tmp_2*tmp_30 - tmp_27*tmp_5);
	out[3] = tmp_1*(-tmp_34*tmp_4 + tmp_35*tmp_7) + tmp_16*tmp_36 + tmp_9*(tmp_11*tmp_34 - tmp_13*tmp_35);
	out[4] = tmp_18*(tmp_20*tmp_34 - tmp_22*tmp_35) + tmp_23*(-tmp_24*tmp_34 + tmp_25*tmp_35) + tmp_26*tmp_36;
	out[5] = tmp_31*tmp_36 + tmp_32*(-tmp_28*tmp_34 + tmp_29*tmp_35) + tmp_33*(-tmp_27*tmp_35 + tmp_30*tmp_34);
	out[6] = tmp_1*(-tmp_37*tmp_4 + tmp_38*tmp_7) + tmp_16*tmp_39 + tmp_9*(tmp_11*tmp_37 - tmp_13*tmp_38);
	out[7] = tmp_18*(tmp_20*tmp_37 - tmp_22*tmp_38) + tmp_23*(-tmp_24*tmp_37 + tmp_25*tmp_38) + tmp_26*tmp_39;
	out[8] = tmp_31*tmp_39 + tmp_32*(-tmp_28*tmp_37 + tmp_29*tmp_38) + tmp_33*(-tmp_27*tmp_38 + tmp_30*tmp_37);
}
}
#undef R
#undef I