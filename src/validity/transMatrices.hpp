#pragma once
#include "numbers/Matrix.hpp"

namespace element_validity {
	template<int n, int s, int p>
	void initMatrices(
		Matrix<Interval> &l2b,
		span<Matrix<Interval>> ssd
	);
	
	template<int n, int s, int p>
	void initMatricesT(
		Matrix<Interval> &l2b,
		span<Matrix<Interval>> tsd,
		span<Matrix<Interval>> ssd
	);
}