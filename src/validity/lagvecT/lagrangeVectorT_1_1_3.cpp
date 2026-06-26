#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVectorT<1, 1, 3>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 8);
	assert(out.size() == 6);
	std::array<Interval, 8> cp;
	for (int i = 0; i < 8; ++i) cp[i] = cpFP[i];
	out[0] = -R(11, 2)*cp[0] + 9*cp[2] - R(9, 2)*cp[4] + cp[6];
	out[1] = -R(11, 2)*cp[1] + 9*cp[3] - R(9, 2)*cp[5] + cp[7];
	out[2] = (R(1, 8))*cp[0] - R(27, 8)*cp[2] + (R(27, 8))*cp[4] - R(1, 8)*cp[6];
	out[3] = (R(1, 8))*cp[1] - R(27, 8)*cp[3] + (R(27, 8))*cp[5] - R(1, 8)*cp[7];
	out[4] = -cp[0] + (R(9, 2))*cp[2] - 9*cp[4] + (R(11, 2))*cp[6];
	out[5] = -cp[1] + (R(9, 2))*cp[3] - 9*cp[5] + (R(11, 2))*cp[7];
}
}
#undef R
#undef I