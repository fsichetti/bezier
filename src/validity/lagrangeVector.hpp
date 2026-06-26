#pragma once
#include "numbers/Interval.hpp"
#include <array>

namespace element_validity {
	template<int n, int s, int p>
	void lagrangeVector(const span<const fp_t>, const span<Interval>);

	template<int n, int s, int p>
	void lagrangeVectorT(const span<const fp_t>, const span<Interval>);
}