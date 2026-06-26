#pragma once
#include <vector>
#include <string>
#include "Interval.hpp"

/* Simple class for precomputed matrices */

namespace element_validity {

template<typename T>
class Matrix{
	private:
	// Store each row as a list of column indices and values 
	using Row = std::vector<std::pair<int, T>>;
	std::vector<Row> store;
	int s;

	private:
	// Write operation that only works in ascending order of column index.
	void pushToRow(int i, int j, const T &v);

	void mult(const span<const T> src, const span<T> dst) const;

	public:
	Matrix() {}
	Matrix(const std::vector<lint> &data) { fill(data); }
	void fill(const std::vector<lint> &data);
	inline const Row &getRow(int i) const { return store.at(i); }
	inline void resize(int size) { s = size; store.resize(s); }

	// Apply matrix to column vector src and directly into dst 
	void mult(const span<const T> src, std::vector<T> &dst) const;

	// Apply matrix to column vector v
	// inline std::vector<T> mult(const std::vector<T> &src) const {
	// 	std::vector<T> res;
	// 	mult(src, res);
	// 	return res;
	// };

	// Stream output
	// template<typename U>
	// friend std::ostream& operator<<(std::ostream& ost, const Matrix<U>& m) {
	// 	const U zero = 0;
	// 	for(int i=0; i<m.s; ++i) {
	// 		const typename Matrix<U>::Row &r = m.getRow(i);
	// 		int j = 0;
	// 		for (const std::pair<int, U> &p : r) {
	// 			while (j < p.first) { ost << zero << '\t'; ++j; }
	// 			ost << p.second << '\t'; ++j;
	// 		}
	// 		while (j < m.s) { ost << zero << '\t'; ++j; }
	// 		ost << std::endl;
	// 	}
	// 	return ost;
	// }
};
}

#include "Matrix.tpp"
