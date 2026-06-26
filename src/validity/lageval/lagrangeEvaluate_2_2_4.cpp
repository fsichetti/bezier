#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<2, 2, 4>(
	const span<const fp_t> xFP,
	const span<const Interval> lagVec
) {
	assert(xFP.size() == 2);
	assert(lagVec.size() == 28);
	std::array<Interval, 2> x;
	for (int i = 0; i < 2; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = x[0]*x[1];
	I tmp_1 = powi(x[0], 2);
	I tmp_2 = powi(x[0], 3);
	I tmp_3 = powi(x[0], 4);
	I tmp_4 = powi(x[0], 5);
	I tmp_5 = powi(x[0], 6);
	I tmp_6 = (R(324, 5))*tmp_5;
	I tmp_7 = powi(x[1], 2);
	I tmp_8 = powi(x[1], 3);
	I tmp_9 = powi(x[1], 4);
	I tmp_10 = powi(x[1], 5);
	I tmp_11 = powi(x[1], 6);
	I tmp_12 = (R(324, 5))*tmp_11;
	I tmp_13 = tmp_7*x[0];
	I tmp_14 = tmp_8*x[0];
	I tmp_15 = tmp_9*x[0];
	I tmp_16 = tmp_10*x[0];
	I tmp_17 = (R(1944, 5))*tmp_16;
	I tmp_18 = tmp_1*x[1];
	I tmp_19 = tmp_2*x[1];
	I tmp_20 = tmp_3*x[1];
	I tmp_21 = tmp_4*x[1];
	I tmp_22 = (R(1944, 5))*tmp_21;
	I tmp_23 = tmp_1*tmp_7;
	I tmp_24 = tmp_1*tmp_8;
	I tmp_25 = tmp_1*tmp_9;
	I tmp_26 = 972*tmp_25;
	I tmp_27 = tmp_2*tmp_7;
	I tmp_28 = tmp_2*tmp_8;
	I tmp_29 = 1296*tmp_28;
	I tmp_30 = tmp_3*tmp_7;
	I tmp_31 = 972*tmp_30;
	I tmp_32 = 1944*tmp_16;
	I tmp_33 = 1944*tmp_30;
	I tmp_34 = (R(1944, 5))*tmp_11;
	I tmp_35 = 3888*tmp_25;
	I tmp_36 = -5184*tmp_2*tmp_7;
	I tmp_37 = 3888*tmp_28;
	I tmp_38 = (R(1566, 5))*tmp_0 + 5022*tmp_23 + tmp_37;
	I tmp_39 = 5832*tmp_25;
	I tmp_40 = 972*tmp_11;
	I tmp_41 = -162*tmp_20 + tmp_31;
	I tmp_42 = 3888*tmp_10*x[0];
	I tmp_43 = tmp_37 + tmp_42;
	I tmp_44 = (R(513, 2))*tmp_0 + 4671*tmp_23;
	I tmp_45 = -5184*tmp_1*tmp_8;
	I tmp_46 = 72*tmp_19 - 648*tmp_27 + tmp_29;
	I tmp_47 = 148*tmp_0 + 1836*tmp_23;
	I tmp_48 = 297*tmp_7;
	I tmp_49 = -27*tmp_18 - 972*tmp_24;
	I tmp_50 = tmp_1*tmp_48;
	I tmp_51 = (R(99, 2))*tmp_0 + tmp_50;
	I tmp_52 = (R(36, 5))*x[0]*x[1];
	I tmp_53 = -648*tmp_15;
	I tmp_54 = -90*tmp_13 + 378*tmp_14 + tmp_17 + tmp_52 + tmp_53;
	I tmp_55 = 1944*tmp_25;
	I tmp_56 = 1944*tmp_21;
	I tmp_57 = 3888*tmp_30;
	I tmp_58 = (R(1944, 5))*tmp_5;
	I tmp_59 = 7776*tmp_25;
	I tmp_60 = 11664*tmp_28;
	I tmp_61 = 7776*tmp_30;
	I tmp_62 = -648*tmp_20;
	I tmp_63 = tmp_57 + tmp_62;
	I tmp_64 = 360*tmp_0 + 11232*tmp_23 + tmp_60;
	I tmp_65 = -1944*tmp_27;
	I tmp_66 = 216*tmp_19 + tmp_65;
	I tmp_67 = 180*tmp_0 + 3996*tmp_23;
	I tmp_68 = -1944*tmp_24;
	I tmp_69 = 54*tmp_0 + 594*tmp_23;
	I tmp_70 = 5832*tmp_30;
	I tmp_71 = 972*tmp_5;
	I tmp_72 = -162*tmp_15 + tmp_26;
	I tmp_73 = 3888*tmp_4*x[1];
	I tmp_74 = tmp_37 + tmp_73;
	I tmp_75 = tmp_35 + tmp_53;
	I tmp_76 = 36*x[1];
	I tmp_77 = 2484*tmp_23 + tmp_37 + tmp_76*x[0];
	I tmp_78 = (R(9, 2))*tmp_0 + tmp_50;
	I tmp_79 = 72*tmp_14 - 648*tmp_24;
	I tmp_80 = 216*tmp_14 + tmp_68;
	I tmp_81 = -27*tmp_13 - 972*tmp_27;
	I tmp_82 = -90*tmp_18 + 378*tmp_19 + tmp_22 + tmp_52 + tmp_62;
	acc += lagVec[0] * ((R(812, 5))*tmp_0 + (R(406, 5))*tmp_1 - R(1134, 5)*tmp_10 + tmp_12 - R(1323, 2)*tmp_13 + 1260*tmp_14 - 1134*tmp_15 + tmp_17 - R(1323, 2)*tmp_18 + 1260*tmp_19 - R(441, 2)*tmp_2 - 1134*tmp_20 + tmp_22 + 1890*tmp_23 - 2268*tmp_24 + tmp_26 - 2268*tmp_27 + tmp_29 + 315*tmp_3 + tmp_31 - R(1134, 5)*tmp_4 + tmp_6 + (R(406, 5))*tmp_7 - R(441, 2)*tmp_8 + 315*tmp_9 - R(147, 10)*x[0] - R(147, 10)*x[1] + Interval(1.));
	acc += lagVec[1] * (7776*tmp_1*tmp_8 + 1044*tmp_1*x[1] + 1296*tmp_10 - 5022*tmp_14 - 1674*tmp_19 - tmp_22 + 1296*tmp_3*x[1] - tmp_32 - tmp_33 - tmp_34 - tmp_35 - tmp_36 - tmp_38 + 2088*tmp_7*x[0] - R(1566, 5)*tmp_7 + 1044*tmp_8 + 5184*tmp_9*x[0] - 1674*tmp_9 + 36*x[1]);
	acc += lagVec[2] * (-3078*tmp_10 - 2610*tmp_13 + 7884*tmp_14 - 9396*tmp_15 - R(1071, 2)*tmp_18 + 486*tmp_19 - 9720*tmp_24 - 3564*tmp_27 + tmp_39 + tmp_40 + tmp_41 + tmp_43 + tmp_44 + (R(1053, 2))*tmp_7 - R(4149, 2)*tmp_8 + 3699*tmp_9 - 45*x[1]);
	acc += lagVec[3] * (180*tmp_1*x[1] + 3888*tmp_10 - 1296*tmp_11 - 6120*tmp_14 - tmp_35 - tmp_42 - tmp_45 - tmp_46 - tmp_47 + 1692*tmp_7*x[0] - 508*tmp_7 + 2232*tmp_8 + 8424*tmp_9*x[0] - 4356*tmp_9 + 40*x[1]);
	acc += lagVec[4] * (-2754*tmp_10 - R(1197, 2)*tmp_13 + 2376*tmp_14 - 3726*tmp_15 + tmp_26 + tmp_32 + tmp_40 + tmp_48 + tmp_49 + tmp_51 - R(2763, 2)*tmp_8 + 2889*tmp_9 - R(45, 2)*x[1]);
	acc += lagVec[5] * ((R(5184, 5))*tmp_10 - tmp_34 - tmp_54 - R(486, 5)*tmp_7 + 468*tmp_8 - 1026*tmp_9 + (R(36, 5))*x[1]);
	acc += lagVec[6] * (-162*tmp_10 + tmp_12 + (R(137, 10))*tmp_7 - R(135, 2)*tmp_8 + 153*tmp_9 - x[1]);
	acc += lagVec[7] * (2088*tmp_1*x[1] - R(1566, 5)*tmp_1 - 1674*tmp_14 - tmp_17 - 5022*tmp_19 + 7776*tmp_2*tmp_7 + 1044*tmp_2 + 5184*tmp_3*x[1] - 1674*tmp_3 - tmp_38 + 1296*tmp_4 - tmp_45 - tmp_55 - tmp_56 - tmp_57 - tmp_58 + 1044*tmp_7*x[0] + 1296*tmp_9*x[0] + 36*x[0]);
	acc += lagVec[8] * (540*tmp_0 - 3078*tmp_13 + 6426*tmp_14 - 5832*tmp_15 - 3078*tmp_18 + 6426*tmp_19 - 5832*tmp_20 + 12852*tmp_23 - 17496*tmp_24 - 17496*tmp_27 + tmp_32 + tmp_56 + tmp_59 + tmp_60 + tmp_61);
	acc += lagVec[9] * (21384*tmp_1*tmp_8 + 1332*tmp_1*x[1] - 9612*tmp_14 - 1620*tmp_19 + 11664*tmp_2*tmp_7 - 11664*tmp_25 - tmp_42 - tmp_63 - tmp_64 + 3492*tmp_7*x[0] + 10368*tmp_9*x[0]);
	acc += lagVec[10] * (-2016*tmp_13 + 7020*tmp_14 - 9072*tmp_15 - 396*tmp_18 - 11016*tmp_24 + tmp_43 + tmp_59 + tmp_66 + tmp_67);
	acc += lagVec[11] * (54*tmp_1*x[1] - 2538*tmp_14 - tmp_32 - tmp_55 - tmp_68 - tmp_69 + 648*tmp_7*x[0] + 3888*tmp_9*x[0]);
	acc += lagVec[12] * (tmp_54);
	acc += lagVec[13] * ((R(1053, 2))*tmp_1 - R(1071, 2)*tmp_13 + 486*tmp_14 - 2610*tmp_18 + 7884*tmp_19 - R(4149, 2)*tmp_2 - 9396*tmp_20 - 3564*tmp_24 - 9720*tmp_27 + 3699*tmp_3 - 3078*tmp_4 + tmp_44 + tmp_70 + tmp_71 + tmp_72 + tmp_74 - 45*x[0]);
	acc += lagVec[14] * (11664*tmp_1*tmp_8 + 3492*tmp_1*x[1] - 1620*tmp_14 - 9612*tmp_19 + 21384*tmp_2*tmp_7 + 10368*tmp_3*x[1] - 11664*tmp_30 - tmp_64 + 1332*tmp_7*x[0] - tmp_73 - tmp_75);
	acc += lagVec[15] * (135*tmp_0 - 1107*tmp_13 + 1944*tmp_14 - 972*tmp_15 - 1107*tmp_18 + 1944*tmp_19 - 972*tmp_20 + 8748*tmp_23 - 13608*tmp_24 - 13608*tmp_27 + tmp_39 + tmp_60 + tmp_70);
	acc += lagVec[16] * (6480*tmp_1*tmp_8 + 252*tmp_1*x[1] - 972*tmp_14 - tmp_66 + 360*tmp_7*x[0] - tmp_75 - tmp_77);
	acc += lagVec[17] * (-R(99, 2)*tmp_13 + 162*tmp_14 + tmp_49 + tmp_72 + tmp_78);
	acc += lagVec[18] * (1692*tmp_1*x[1] - 508*tmp_1 - 6120*tmp_19 + 2232*tmp_2 - tmp_29 + 8424*tmp_3*x[1] - 4356*tmp_3 - tmp_36 + 3888*tmp_4 - tmp_47 - 1296*tmp_5 - tmp_57 + 180*tmp_7*x[0] - tmp_73 - tmp_79 + 40*x[0]);
	acc += lagVec[19] * (-396*tmp_13 - 2016*tmp_18 + 7020*tmp_19 - 9072*tmp_20 - 11016*tmp_27 + tmp_61 + tmp_67 + tmp_74 + tmp_80);
	acc += lagVec[20] * (360*tmp_1*x[1] - 972*tmp_19 + 6480*tmp_2*tmp_7 - tmp_63 + 252*tmp_7*x[0] - tmp_77 - tmp_80);
	acc += lagVec[21] * (4*tmp_0 - tmp_1*tmp_76 + 324*tmp_23 + tmp_46 - 36*tmp_7*x[0] + tmp_79);
	acc += lagVec[22] * (297*tmp_1 - R(1197, 2)*tmp_18 + 2376*tmp_19 - R(2763, 2)*tmp_2 - 3726*tmp_20 + 2889*tmp_3 + tmp_31 - 2754*tmp_4 + tmp_51 + tmp_56 + tmp_71 + tmp_81 - R(45, 2)*x[0]);
	acc += lagVec[23] * (648*tmp_1*x[1] - 2538*tmp_19 + 3888*tmp_3*x[1] - tmp_33 - tmp_56 - tmp_65 - tmp_69 + 54*tmp_7*x[0]);
	acc += lagVec[24] * (-R(99, 2)*tmp_18 + 162*tmp_19 + tmp_41 + tmp_78 + tmp_81);
	acc += lagVec[25] * (-R(486, 5)*tmp_1 + 468*tmp_2 - 1026*tmp_3 + (R(5184, 5))*tmp_4 - tmp_58 - tmp_82 + (R(36, 5))*x[0]);
	acc += lagVec[26] * (tmp_82);
	acc += lagVec[27] * ((R(137, 10))*tmp_1 - R(135, 2)*tmp_2 + 153*tmp_3 - 162*tmp_4 + tmp_6 - x[0]);
	return acc;
}}
#undef R
#undef I