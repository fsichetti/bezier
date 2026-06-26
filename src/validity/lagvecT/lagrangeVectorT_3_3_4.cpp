#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

#include "lagrangeVectorT_3_3_4_chunks.hpp"

namespace element_validity {
template<>
void lagrangeVectorT<3, 3, 4>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 210);
	assert(out.size() == 880);
	std::array<Interval, 210> cp;
	for (int i = 0; i < 210; ++i) cp[i] = cpFP[i];
	chunkT_3_3_4_0(cp, out);
	chunkT_3_3_4_1(cp, out);
	chunkT_3_3_4_2(cp, out);
	chunkT_3_3_4_3(cp, out);
	chunkT_3_3_4_4(cp, out);
	chunkT_3_3_4_5(cp, out);
	chunkT_3_3_4_6(cp, out);
	chunkT_3_3_4_7(cp, out);
	chunkT_3_3_4_8(cp, out);
	chunkT_3_3_4_9(cp, out);
	chunkT_3_3_4_10(cp, out);
	chunkT_3_3_4_11(cp, out);
	chunkT_3_3_4_12(cp, out);
	chunkT_3_3_4_13(cp, out);
	chunkT_3_3_4_14(cp, out);
	chunkT_3_3_4_15(cp, out);
	chunkT_3_3_4_16(cp, out);
	chunkT_3_3_4_17(cp, out);
}
}
#undef R
#undef I