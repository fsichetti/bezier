#pragma once
#include "utils/globals.hpp"
#include <vector>

namespace element_validity {
	template<int n, int s, int p>
	void cornerIndices(std::vector<int> &);

	template<int n, int s, int p>
	void cornerIndicesT(std::vector<int> &);
}