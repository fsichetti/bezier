#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVectorT<2, 1, 1>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 16);
	assert(out.size() == 12);
	std::array<Interval, 16> cp;
	for (int i = 0; i < 16; ++i) cp[i] = cpFP[i];
	I tmp_0 = -cp[0] + cp[12];
	I tmp_1 = -cp[6];
	I tmp_2 = -cp[2] - tmp_1;
	I tmp_3 = -cp[0] + cp[4];
	I tmp_4 = cp[14] - cp[2];
	I tmp_5 = (R(1, 2))*cp[12];
	I tmp_6 = (R(1, 2))*cp[13];
	I tmp_7 = (R(1, 2))*cp[0] + (R(1, 2))*cp[1];
	I tmp_8 = tmp_5 + tmp_6 - tmp_7;
	I tmp_9 = (R(1, 2))*cp[2];
	I tmp_10 = (R(1, 2))*cp[3];
	I tmp_11 = -R(1, 2)*cp[6] - R(1, 2)*cp[7];
	I tmp_12 = -tmp_10 - tmp_11 - tmp_9;
	I tmp_13 = (R(1, 2))*cp[4];
	I tmp_14 = (R(1, 2))*cp[5];
	I tmp_15 = tmp_13 + tmp_14 - tmp_7;
	I tmp_16 = (R(1, 2))*cp[14];
	I tmp_17 = (R(1, 2))*cp[15];
	I tmp_18 = -tmp_10 + tmp_16 + tmp_17 - tmp_9;
	I tmp_19 = cp[13] - cp[1];
	I tmp_20 = -cp[7];
	I tmp_21 = -cp[3] - tmp_20;
	I tmp_22 = cp[15] - cp[3];
	I tmp_23 = -cp[1] + cp[5];
	I tmp_24 = cp[10] - cp[14];
	I tmp_25 = -cp[8];
	I tmp_26 = -cp[12] - tmp_25;
	I tmp_27 = (R(1, 2))*cp[10] + (R(1, 2))*cp[11];
	I tmp_28 = -tmp_16 - tmp_17 + tmp_27;
	I tmp_29 = -R(1, 2)*cp[8] - R(1, 2)*cp[9];
	I tmp_30 = -tmp_29 - tmp_5 - tmp_6;
	I tmp_31 = cp[11] - cp[15];
	I tmp_32 = -cp[9];
	I tmp_33 = -cp[13] - tmp_32;
	I tmp_34 = cp[10] + tmp_1;
	I tmp_35 = -cp[4] - tmp_25;
	I tmp_36 = tmp_11 + tmp_27;
	I tmp_37 = -tmp_13 - tmp_14 - tmp_29;
	I tmp_38 = cp[11] + tmp_20;
	I tmp_39 = -cp[5] - tmp_32;
	out[0] = -tmp_0*tmp_2 + tmp_3*tmp_4;
	out[1] = -tmp_12*tmp_8 + tmp_15*tmp_18;
	out[2] = -tmp_19*tmp_21 + tmp_22*tmp_23;
	out[3] = -tmp_0*tmp_24 + tmp_26*tmp_4;
	out[4] = tmp_18*tmp_30 - tmp_28*tmp_8;
	out[5] = -tmp_19*tmp_31 + tmp_22*tmp_33;
	out[6] = -tmp_2*tmp_35 + tmp_3*tmp_34;
	out[7] = -tmp_12*tmp_37 + tmp_15*tmp_36;
	out[8] = -tmp_21*tmp_39 + tmp_23*tmp_38;
	out[9] = -tmp_24*tmp_35 + tmp_26*tmp_34;
	out[10] = -tmp_28*tmp_37 + tmp_30*tmp_36;
	out[11] = -tmp_31*tmp_39 + tmp_33*tmp_38;
}
}
#undef R
#undef I