#pragma once
#include "utils/globals.hpp"

namespace element_validity {

	// Binomial coefficient
	constexpr int binom(int n, int k) {
		// edge cases
		if (k>n) return 0;
		if (k==0 || k==n) return 1;
		k = std::min(k, n-k);
		// accumulate results
		int acc = 1;
		for (int i=0; i<k; ++i) {
			acc = (acc * (n-i)) / (i+1);
		}
		return acc;
	}

	// Positive integer power
	template<typename T>
	constexpr T powi(T n, int e) {
		switch(e) {
		// Base cases
		case 0: return 1;
		case 1: return n;
		default:
			// Exponent is even: divide by 2 and recurse
			if ((e % 2) == 0) {
				const T a = powi(n, e >> 1);
				return a*a;
			}
			// Exponent is not divisible by 2
			return n * powi(n, e - 1);
		}
	}

	// Number of control points
	constexpr int nControlGeoMap(int n, int s, int p) {
		return binom(p+s, s) * powi(p+1, n - s);
	}

	// Number of control points for the jacobian
	constexpr int nControlJacobian(int n, int s, int p, bool t) {
		int L = binom(n*p, s) * powi(n * p, n - s);
		return t ? L * (n+1) : L;
	}
}