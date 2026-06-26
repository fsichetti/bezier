#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<1, 1, 4>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 5);
	assert(out.size() == 4);
	std::array<Interval, 5> cp;
	for (int i = 0; i < 5; ++i) cp[i] = cpFP[i];
	I tmp_0 = 12*cp[2];
	I tmp_1 = -R(148, 27)*cp[2];
	out[0] = -R(25, 3)*cp[0] + 16*cp[1] + (R(16, 3))*cp[3] - cp[4] - tmp_0;
	out[1] = -R(7, 81)*cp[0] - R(368, 81)*cp[1] - R(80, 81)*cp[3] + (R(11, 81))*cp[4] - tmp_1;
	out[2] = -R(11, 81)*cp[0] + (R(80, 81))*cp[1] + (R(368, 81))*cp[3] + (R(7, 81))*cp[4] + tmp_1;
	out[3] = cp[0] - R(16, 3)*cp[1] - 16*cp[3] + (R(25, 3))*cp[4] + tmp_0;
}
}
#undef R
#undef I