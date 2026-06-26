#ifdef EIGEN_INTERFACE
#include "eigen.hpp"
#include <cassert>

namespace element_validity {
std::vector<double> convertEigenMatrix(const Eigen::MatrixXd& mat) {
	std::vector<double> vec;
	vec.reserve(mat.size());
	for (int i = 0; i < mat.rows(); ++i) {
		for (int j = 0; j < mat.cols(); ++j) {
			vec.push_back(mat(i, j));
		}
	}
	return vec;
}
std::vector<double> convertEigenMatrix(
	const Eigen::MatrixXd& mat1,
	const Eigen::MatrixXd& mat2
) {
	std::vector<double> vec;
	assert(mat1.rows() == mat2.rows() && mat1.cols() == mat2.cols());
	vec.reserve(mat1.size() * 2);
	for (int i = 0; i < mat1.rows(); ++i) {
		for (int j = 0; j < mat1.cols(); ++j) {
			vec.push_back(mat1(i, j));
			vec.push_back(mat2(i, j));
		}
	}
	return vec;
}
}
#endif