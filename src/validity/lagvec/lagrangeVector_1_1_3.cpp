#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<1, 1, 3>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 4);
	assert(out.size() == 3);
	std::array<Interval, 4> cp;
	for (int i = 0; i < 4; ++i) cp[i] = cpFP[i];
	out[0] = -R(11, 2)*cp[0] + 9*cp[1] - R(9, 2)*cp[2] + cp[3];
	out[1] = (R(1, 8))*cp[0] - R(27, 8)*cp[1] + (R(27, 8))*cp[2] - R(1, 8)*cp[3];
	out[2] = -cp[0] + (R(9, 2))*cp[1] - 9*cp[2] + (R(11, 2))*cp[3];
}
}
#undef R
#undef I