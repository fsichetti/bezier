#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
void lagrangeVectorT<1, 1, 5>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 12);
	assert(out.size() == 10);
	std::array<Interval, 12> cp;
	for (int i = 0; i < 12; ++i) cp[i] = cpFP[i];
	out[0] = -R(137, 12)*cp[0] + cp[10] + 25*cp[2] - 25*cp[4] + (R(50, 3))*cp[6] - R(25, 4)*cp[8];
	out[1] = cp[11] - R(137, 12)*cp[1] + 25*cp[3] - 25*cp[5] + (R(50, 3))*cp[7] - R(25, 4)*cp[9];
	out[2] = -R(423, 2048)*cp[0] - R(731, 6144)*cp[10] - R(12525, 2048)*cp[2] + (R(25975, 3072))*cp[4] - R(2925, 1024)*cp[6] + (R(1725, 2048))*cp[8];
	out[3] = -R(731, 6144)*cp[11] - R(423, 2048)*cp[1] - R(12525, 2048)*cp[3] + (R(25975, 3072))*cp[5] - R(2925, 1024)*cp[7] + (R(1725, 2048))*cp[9];
	out[4] = -R(3, 128)*cp[0] + (R(3, 128))*cp[10] + (R(125, 384))*cp[2] - R(375, 64)*cp[4] + (R(375, 64))*cp[6] - R(125, 384)*cp[8];
	out[5] = (R(3, 128))*cp[11] - R(3, 128)*cp[1] + (R(125, 384))*cp[3] - R(375, 64)*cp[5] + (R(375, 64))*cp[7] - R(125, 384)*cp[9];
	out[6] = (R(731, 6144))*cp[0] + (R(423, 2048))*cp[10] - R(1725, 2048)*cp[2] + (R(2925, 1024))*cp[4] - R(25975, 3072)*cp[6] + (R(12525, 2048))*cp[8];
	out[7] = (R(423, 2048))*cp[11] + (R(731, 6144))*cp[1] - R(1725, 2048)*cp[3] + (R(2925, 1024))*cp[5] - R(25975, 3072)*cp[7] + (R(12525, 2048))*cp[9];
	out[8] = -cp[0] + (R(137, 12))*cp[10] + (R(25, 4))*cp[2] - R(50, 3)*cp[4] + 25*cp[6] - 25*cp[8];
	out[9] = (R(137, 12))*cp[11] - cp[1] + (R(25, 4))*cp[3] - R(50, 3)*cp[5] + 25*cp[7] - 25*cp[9];
}
}
#undef R
#undef I