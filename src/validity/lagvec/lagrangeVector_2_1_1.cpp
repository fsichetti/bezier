#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<2, 1, 1>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 8);
	assert(out.size() == 4);
	std::array<Interval, 8> cp;
	for (int i = 0; i < 8; ++i) cp[i] = cpFP[i];
	I tmp_0 = -cp[0] + cp[2];
	I tmp_1 = -cp[7];
	I tmp_2 = -cp[1] - tmp_1;
	I tmp_3 = -cp[6];
	I tmp_4 = -cp[0] - tmp_3;
	I tmp_5 = -cp[1] + cp[3];
	I tmp_6 = cp[5] + tmp_1;
	I tmp_7 = cp[4] + tmp_3;
	I tmp_8 = -cp[3] + cp[5];
	I tmp_9 = -cp[2] + cp[4];
	out[0] = tmp_0*tmp_2 - tmp_4*tmp_5;
	out[1] = tmp_2*tmp_7 - tmp_4*tmp_6;
	out[2] = tmp_0*tmp_8 - tmp_5*tmp_9;
	out[3] = -tmp_6*tmp_9 + tmp_7*tmp_8;
}
}
#undef R
#undef I