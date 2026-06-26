#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVector<1, 1, 5>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 6);
	assert(out.size() == 5);
	std::array<Interval, 6> cp;
	for (int i = 0; i < 6; ++i) cp[i] = cpFP[i];
	out[0] = -R(137, 12)*cp[0] + 25*cp[1] - 25*cp[2] + (R(50, 3))*cp[3] - R(25, 4)*cp[4] + cp[5];
	out[1] = -R(423, 2048)*cp[0] - R(12525, 2048)*cp[1] + (R(25975, 3072))*cp[2] - R(2925, 1024)*cp[3] + (R(1725, 2048))*cp[4] - R(731, 6144)*cp[5];
	out[2] = -R(3, 128)*cp[0] + (R(125, 384))*cp[1] - R(375, 64)*cp[2] + (R(375, 64))*cp[3] - R(125, 384)*cp[4] + (R(3, 128))*cp[5];
	out[3] = (R(731, 6144))*cp[0] - R(1725, 2048)*cp[1] + (R(2925, 1024))*cp[2] - R(25975, 3072)*cp[3] + (R(12525, 2048))*cp[4] + (R(423, 2048))*cp[5];
	out[4] = -cp[0] + (R(25, 4))*cp[1] - R(50, 3)*cp[2] + 25*cp[3] - 25*cp[4] + (R(137, 12))*cp[5];
}
}
#undef R
#undef I