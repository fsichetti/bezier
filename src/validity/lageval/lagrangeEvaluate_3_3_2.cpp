#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<3, 3, 2>(
	const span<const fp_t> xFP,
	const span<const Interval> lagVec
) {
	assert(xFP.size() == 3);
	assert(lagVec.size() == 20);
	std::array<Interval, 3> x;
	for (int i = 0; i < 3; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = powi(x[2], 3);
	I tmp_1 = (R(9, 2))*tmp_0;
	I tmp_2 = powi(x[1], 3);
	I tmp_3 = (R(9, 2))*tmp_2;
	I tmp_4 = powi(x[0], 3);
	I tmp_5 = (R(9, 2))*tmp_4;
	I tmp_6 = x[1]*x[2];
	I tmp_7 = 27*tmp_6*x[0];
	I tmp_8 = powi(x[0], 2);
	I tmp_9 = powi(x[1], 2);
	I tmp_10 = powi(x[2], 2);
	I tmp_11 = (R(27, 2))*x[2];
	I tmp_12 = tmp_11*tmp_9;
	I tmp_13 = tmp_11*tmp_8;
	I tmp_14 = tmp_12 + tmp_13;
	I tmp_15 = (R(27, 2))*x[1];
	I tmp_16 = tmp_10*tmp_15;
	I tmp_17 = tmp_15*tmp_8;
	I tmp_18 = tmp_16 + tmp_17;
	I tmp_19 = (R(27, 2))*x[0];
	I tmp_20 = tmp_10*tmp_19;
	I tmp_21 = tmp_19*tmp_9;
	I tmp_22 = tmp_20 + tmp_21;
	I tmp_23 = 27*tmp_10;
	I tmp_24 = tmp_23*x[0];
	I tmp_25 = (R(27, 2))*tmp_0;
	I tmp_26 = (R(45, 2))*x[0];
	I tmp_27 = -tmp_26*x[2];
	I tmp_28 = -R(45, 2)*tmp_6;
	I tmp_29 = tmp_23*x[1] + tmp_7;
	I tmp_30 = (R(9, 2))*x[2];
	I tmp_31 = -tmp_30*x[1];
	I tmp_32 = tmp_16 + tmp_31;
	I tmp_33 = -tmp_30*x[0];
	I tmp_34 = tmp_20 + tmp_33;
	I tmp_35 = 27*tmp_9;
	I tmp_36 = tmp_35*x[2];
	I tmp_37 = (R(27, 2))*tmp_2;
	I tmp_38 = -tmp_26*x[1];
	I tmp_39 = tmp_35*x[0] + tmp_7;
	I tmp_40 = (R(9, 2))*x[1];
	I tmp_41 = tmp_12 + tmp_31;
	I tmp_42 = -tmp_40*x[0];
	I tmp_43 = tmp_21 + tmp_42;
	I tmp_44 = 27*tmp_8;
	I tmp_45 = tmp_44*x[1];
	I tmp_46 = (R(27, 2))*tmp_4;
	I tmp_47 = tmp_44*x[2] + tmp_7;
	I tmp_48 = tmp_13 + tmp_33;
	I tmp_49 = tmp_17 + tmp_42;
	acc += lagVec[0] * (-tmp_1 + 9*tmp_10 - tmp_14 - tmp_18 - tmp_22 - tmp_3 - tmp_5 - tmp_7 + 9*tmp_8 + 9*tmp_9 + 18*x[0]*x[1] + 18*x[0]*x[2] - R(11, 2)*x[0] + 18*x[1]*x[2] - R(11, 2)*x[1] - R(11, 2)*x[2] + Interval(1));
	acc += lagVec[1] * (-R(45, 2)*tmp_10 + tmp_14 + tmp_24 + tmp_25 + tmp_27 + tmp_28 + tmp_29 + 9*x[2]);
	acc += lagVec[2] * (18*tmp_10 - tmp_25 - tmp_30 - tmp_32 - tmp_34);
	acc += lagVec[3] * (tmp_1 - R(9, 2)*tmp_10 + x[2]);
	acc += lagVec[4] * (tmp_18 + tmp_28 + tmp_36 + tmp_37 + tmp_38 + tmp_39 - R(45, 2)*tmp_9 + 9*x[1]);
	acc += lagVec[5] * (-tmp_29 - tmp_36 + 27*x[1]*x[2]);
	acc += lagVec[6] * (tmp_32);
	acc += lagVec[7] * (-tmp_37 - tmp_40 - tmp_41 - tmp_43 + 18*tmp_9);
	acc += lagVec[8] * (tmp_41);
	acc += lagVec[9] * (tmp_3 - R(9, 2)*tmp_9 + x[1]);
	acc += lagVec[10] * (tmp_22 + tmp_27 + tmp_38 + tmp_45 + tmp_46 + tmp_47 - R(45, 2)*tmp_8 + 9*x[0]);
	acc += lagVec[11] * (-tmp_24 - tmp_47 + 27*x[0]*x[2]);
	acc += lagVec[12] * (tmp_34);
	acc += lagVec[13] * (-tmp_39 - tmp_45 + 27*x[0]*x[1]);
	acc += lagVec[14] * (tmp_7);
	acc += lagVec[15] * (tmp_43);
	acc += lagVec[16] * (-tmp_46 - tmp_48 - tmp_49 + 18*tmp_8 - R(9, 2)*x[0]);
	acc += lagVec[17] * (tmp_48);
	acc += lagVec[18] * (tmp_49);
	acc += lagVec[19] * (tmp_5 - R(9, 2)*tmp_8 + x[0]);
	return acc;
}}
#undef R
#undef I