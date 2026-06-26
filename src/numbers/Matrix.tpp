#pragma once

#include "Matrix.hpp"

namespace element_validity {

template<typename T>
void Matrix<T>::pushToRow(int i, int j, const T &v) {
	if (store.at(i).size() == 0 || store.at(i).back().first < j)
		store.at(i).emplace_back(j, v);
	else throw std::runtime_error("Attempted out of order insert");
}

/*
Matrices are read from vectors. The first entry is the matrix size,
followed by sequences of {i, j, numerator, denominator} as integers.
*/
template<typename T>
void Matrix<T>::fill(const std::vector<lint> &data) {
	int s = data.size();
	assert(s>0);
	resize(data.at(0));
	for(int k=1; k<s; k+=4) {
		assert(data.at(k+0) >= 0 && data.at(k+1) >= 0);
		pushToRow(
			data.at(k+0),
			data.at(k+1),
			static_cast<T>(data.at(k+2)) / data.at(k+3)
		);
	}
}

template<typename T>
void Matrix<T>::mult(const span<const T> src, const span<T> dst) const {
	assert(src.size() == s);
	assert(dst.size() == s);
	for(int i=0; i<s; ++i) {
		T acc = 0;
		for(const auto &p : getRow(i)) acc += p.second * src[p.first];
		dst[i] = acc;
	}
}

template<typename T>
void Matrix<T>::mult(const span<const T> src, std::vector<T> &dst) const {
	dst.resize(s);
	// mult(src, span<T>(dst));
	assert(src.size() == s);
	assert(dst.size() == s);
	for(int i=0; i<s; ++i) {
		T acc = 0;
		for(const auto &p : getRow(i)) acc += p.second * src[p.first];
		dst[i] = acc;
	}
}
}