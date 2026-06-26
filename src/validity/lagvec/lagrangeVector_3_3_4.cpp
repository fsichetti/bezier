#include "../lagrangeVector.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

#include "lagrangeVector_3_3_4_chunks.hpp"

namespace element_validity {
template<>
void lagrangeVector<3, 3, 4>(const span<const fp_t> cpFP, const span<Interval> out) {
	assert(cpFP.size() == 105);
	assert(out.size() == 220);
	std::array<Interval, 105> cp;
	for (int i = 0; i < 105; ++i) cp[i] = cpFP[i];
	chunk_3_3_4_0(cp, out);
	chunk_3_3_4_1(cp, out);
	chunk_3_3_4_2(cp, out);
	chunk_3_3_4_3(cp, out);
	chunk_3_3_4_4(cp, out);
}
}
#undef R
#undef I