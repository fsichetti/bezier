#include "cornerIndices.hpp"

namespace element_validity {
template<>
void cornerIndices<1, 1, 1>(std::vector<int> &v) { v = {0}; }

template<>
void cornerIndicesT<1, 1, 1>(std::vector<int> &v) { v = {1}; }

template<>
void cornerIndices<1, 1, 2>(std::vector<int> &v) { v = {0,1}; }

template<>
void cornerIndicesT<1, 1, 2>(std::vector<int> &v) { v = {1,3}; }

template<>
void cornerIndices<1, 1, 3>(std::vector<int> &v) { v = {0,2}; }

template<>
void cornerIndicesT<1, 1, 3>(std::vector<int> &v) { v = {1,5}; }

template<>
void cornerIndices<1, 1, 4>(std::vector<int> &v) { v = {0,3}; }

template<>
void cornerIndicesT<1, 1, 4>(std::vector<int> &v) { v = {1,7}; }

template<>
void cornerIndices<1, 1, 5>(std::vector<int> &v) { v = {0,4}; }

template<>
void cornerIndicesT<1, 1, 5>(std::vector<int> &v) { v = {1,9}; }

template<>
void cornerIndices<2, 1, 1>(std::vector<int> &v) { v = {0,1,2,3}; }

template<>
void cornerIndicesT<2, 1, 1>(std::vector<int> &v) { v = {2,5,8,11}; }

template<>
void cornerIndices<2, 1, 2>(std::vector<int> &v) { v = {0,3,12,15}; }

template<>
void cornerIndicesT<2, 1, 2>(std::vector<int> &v) { v = {2,11,38,47}; }

template<>
void cornerIndices<2, 2, 1>(std::vector<int> &v) { v = {0}; }

template<>
void cornerIndicesT<2, 2, 1>(std::vector<int> &v) { v = {2}; }

template<>
void cornerIndices<2, 2, 2>(std::vector<int> &v) { v = {0,2,5}; }

template<>
void cornerIndicesT<2, 2, 2>(std::vector<int> &v) { v = {2,8,17}; }

template<>
void cornerIndices<2, 2, 3>(std::vector<int> &v) { v = {0,4,14}; }

template<>
void cornerIndicesT<2, 2, 3>(std::vector<int> &v) { v = {2,14,44}; }

template<>
void cornerIndices<3, 1, 1>(std::vector<int> &v) { v = {0,2,6,8,18,20,24,26}; }

template<>
void cornerIndicesT<3, 1, 1>(std::vector<int> &v) { v = {3,11,27,35,75,83,99,107}; }

template<>
void cornerIndices<3, 2, 1>(std::vector<int> &v) { v = {0,2,3,5,6,8}; }

template<>
void cornerIndicesT<3, 2, 1>(std::vector<int> &v) { v = {3,11,15,23,27,35}; }

template<>
void cornerIndices<3, 3, 1>(std::vector<int> &v) { v = {0}; }

template<>
void cornerIndicesT<3, 3, 1>(std::vector<int> &v) { v = {3}; }

template<>
void cornerIndices<3, 3, 2>(std::vector<int> &v) { v = {0,3,9,19}; }

template<>
void cornerIndicesT<3, 3, 2>(std::vector<int> &v) { v = {3,15,39,79}; }

template<>
void cornerIndices<2, 2, 4>(std::vector<int> &v) { v = {0,6,27}; }

template<>
void cornerIndicesT<2, 2, 4>(std::vector<int> &v) { v = {2,20,83}; }

template<>
void cornerIndices<3, 3, 3>(std::vector<int> &v) { v = {0,6,27,83}; }

template<>
void cornerIndicesT<3, 3, 3>(std::vector<int> &v) { v = {3,27,111,335}; }

#ifdef INSTANTIATE_3_3_4
template<>
void cornerIndices<3, 3, 4>(std::vector<int> &v) { v = {0,9,54,219}; }

template<>
void cornerIndicesT<3, 3, 4>(std::vector<int> &v) { v = {3,39,219,879}; }
#endif

}
