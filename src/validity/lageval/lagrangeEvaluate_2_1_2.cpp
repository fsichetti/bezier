#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<2, 1, 2>(
	const span<const fp_t> xFP,
	const span<const Interval> lagVec
) {
	assert(xFP.size() == 2);
	assert(lagVec.size() == 16);
	std::array<Interval, 2> x;
	for (int i = 0; i < 2; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = powi(x[1], 2);
	I tmp_1 = powi(x[0], 2);
	I tmp_2 = 81*tmp_1;
	I tmp_3 = tmp_0*tmp_2;
	I tmp_4 = (R(11, 2))*x[0];
	I tmp_5 = tmp_0*x[0];
	I tmp_6 = tmp_1*x[1];
	I tmp_7 = x[0]*x[1];
	I tmp_8 = powi(x[1], 3);
	I tmp_9 = (R(9, 2))*tmp_8;
	I tmp_10 = tmp_1*tmp_8;
	I tmp_11 = tmp_8*x[0];
	I tmp_12 = -R(81, 2)*tmp_10 + (R(99, 4))*tmp_11 - tmp_9;
	I tmp_13 = powi(x[0], 3);
	I tmp_14 = tmp_13*tmp_8;
	I tmp_15 = (R(81, 4))*tmp_14;
	I tmp_16 = (R(9, 2))*tmp_13;
	I tmp_17 = (R(81, 2))*tmp_13;
	I tmp_18 = tmp_13*x[1];
	I tmp_19 = -tmp_0*tmp_17 + tmp_15 - tmp_16 + (R(99, 4))*tmp_18;
	I tmp_20 = -tmp_2*x[1];
	I tmp_21 = (R(45, 2))*tmp_0;
	I tmp_22 = tmp_0*tmp_13;
	I tmp_23 = tmp_17*x[1] - R(405, 4)*tmp_22;
	I tmp_24 = (R(27, 2))*tmp_8;
	I tmp_25 = -R(243, 2)*tmp_10 + (R(297, 4))*tmp_11 - tmp_24;
	I tmp_26 = (R(243, 4))*tmp_14;
	I tmp_27 = tmp_0*tmp_1;
	I tmp_28 = tmp_26 + (R(405, 2))*tmp_27 + (R(99, 2))*tmp_7;
	I tmp_29 = (R(9, 2))*x[1];
	I tmp_30 = tmp_26 - R(81, 2)*tmp_6;
	I tmp_31 = (R(81, 4))*tmp_18 - 81*tmp_22;
	I tmp_32 = 162*tmp_27 + (R(99, 4))*tmp_7;
	I tmp_33 = (R(9, 2))*tmp_0;
	I tmp_34 = tmp_15 + tmp_16*x[1] - R(81, 4)*tmp_22;
	I tmp_35 = (R(81, 2))*tmp_27 + tmp_4*x[1];
	I tmp_36 = -81*tmp_5;
	I tmp_37 = (R(45, 2))*tmp_1;
	I tmp_38 = -R(405, 4)*tmp_10 + (R(81, 2))*tmp_11;
	I tmp_39 = (R(27, 2))*tmp_13;
	I tmp_40 = (R(297, 4))*tmp_18 - R(243, 2)*tmp_22 - tmp_39;
	I tmp_41 = (R(729, 4))*tmp_14;
	I tmp_42 = -R(1215, 4)*tmp_10 + (R(243, 2))*tmp_11 + tmp_41;
	I tmp_43 = (R(243, 2))*tmp_18 - R(1215, 4)*tmp_22;
	I tmp_44 = 405*tmp_27 + (R(81, 2))*tmp_7;
	I tmp_45 = (R(243, 4))*tmp_18 - 243*tmp_22;
	I tmp_46 = tmp_26 - R(81, 2)*tmp_5;
	I tmp_47 = -R(243, 4)*tmp_22 + tmp_39*x[1];
	I tmp_48 = (R(405, 4))*tmp_27 + 9*x[0]*x[1];
	I tmp_49 = (R(81, 4))*tmp_11 - tmp_2*tmp_8;
	I tmp_50 = -243*tmp_10 + (R(243, 4))*tmp_11 + tmp_41;
	I tmp_51 = tmp_26 + tmp_29*x[0] + tmp_3;
	I tmp_52 = -R(81, 4)*tmp_10 + tmp_9*x[0];
	I tmp_53 = -R(243, 4)*tmp_10 + tmp_24*x[0];
	acc += lagVec[0] * (9*tmp_0 + 9*tmp_1 + tmp_12 + tmp_19 + tmp_3 - tmp_4 - R(99, 2)*tmp_5 - R(99, 2)*tmp_6 + (R(121, 4))*tmp_7 - R(11, 2)*x[1] + 1);
	acc += lagVec[1] * ((R(495, 4))*tmp_0*x[0] - tmp_20 - tmp_21 - tmp_23 - tmp_25 - tmp_28 + 9*x[1]);
	acc += lagVec[2] * (18*tmp_0 + tmp_25 - tmp_29 + tmp_30 + tmp_31 + tmp_32 - 99*tmp_5);
	acc += lagVec[3] * ((R(99, 4))*tmp_0*x[0] + 9*tmp_1*x[1] - tmp_12 - tmp_33 - tmp_34 - tmp_35 + x[1]);
	acc += lagVec[4] * ((R(495, 4))*tmp_1*x[1] - tmp_28 - tmp_36 - tmp_37 - tmp_38 - tmp_40 + 9*x[0]);
	acc += lagVec[5] * ((R(2025, 4))*tmp_27 + tmp_42 + tmp_43 - R(405, 2)*tmp_5 - R(405, 2)*tmp_6 + 81*tmp_7);
	acc += lagVec[6] * (162*tmp_0*x[0] + (R(405, 4))*tmp_1*x[1] - tmp_42 - tmp_44 - tmp_45);
	acc += lagVec[7] * (-tmp_37*x[1] + tmp_38 + tmp_46 + tmp_47 + tmp_48);
	acc += lagVec[8] * (18*tmp_1 + tmp_32 + tmp_40 + tmp_46 + tmp_49 - 99*tmp_6 - R(9, 2)*x[0]);
	acc += lagVec[9] * ((R(405, 4))*tmp_0*x[0] + 162*tmp_1*x[1] - tmp_43 - tmp_44 - tmp_50);
	acc += lagVec[10] * (tmp_20 + 324*tmp_27 + tmp_36 + tmp_45 + tmp_50 + (R(81, 4))*tmp_7);
	acc += lagVec[11] * ((R(81, 4))*tmp_0*x[0] + 18*tmp_1*x[1] - tmp_47 - tmp_49 - tmp_51);
	acc += lagVec[12] * (9*tmp_0*x[0] + (R(99, 4))*tmp_1*x[1] - R(9, 2)*tmp_1 - tmp_19 - tmp_35 - tmp_52 + x[0]);
	acc += lagVec[13] * (-tmp_21*x[0] + tmp_23 + tmp_30 + tmp_48 + tmp_53);
	acc += lagVec[14] * (18*tmp_0*x[0] + (R(81, 4))*tmp_1*x[1] - tmp_31 - tmp_51 - tmp_53);
	acc += lagVec[15] * (-tmp_1*tmp_29 + (R(81, 4))*tmp_27 - tmp_33*x[0] + tmp_34 + tmp_52 + tmp_7);
	return acc;
}}
#undef R
#undef I