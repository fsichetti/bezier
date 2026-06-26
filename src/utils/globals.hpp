#pragma once
#include <limits>
#include <iomanip>
#include <cassert>
#include <iostream>
#if __cplusplus >= 202002L
#include <span>
#else
#include "cxx17span.hpp"
#endif

/*
Things that should be available everywhere in the project
*/

namespace element_validity {
	// long integer
	using lint = long;

	// floating point type
	using fp_t = double;
	const auto fp_fmt =
		std::setprecision(std::numeric_limits<fp_t>::max_digits10);

	// validity values
	enum class Validity {valid, invalid, uncertain};

	// span
	template <typename T>
	#if __cplusplus >= 202002L
	using span = std::span<T>;
	#else
	using span = cxx17span<T>;
	#endif
}
