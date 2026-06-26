#pragma once
#include "utils/combinatorics.hpp"
#include "utils/Timer.hpp"
#include "utils/eigen.hpp"
#include "transMatrices.hpp"
#include "lagrangeVector.hpp"
#include "cornerIndices.hpp"

#include <queue>
#include <memory>

namespace element_validity {
class Validator {
	protected:
	Matrix<Interval> matL2B;
	std::vector<int> interpIndices;

	const int nThreads;

	fp_t epsilon = 0.;
	int maxSubdiv = 0;

	struct Subdomain;
	friend std::ostream& operator<<(std::ostream& ost, const Subdomain &s);

	Validator(int nThreads = 1) : nThreads(nThreads) {}

	RealInterval minclusion(const std::vector<Interval> &B) const;

	public:
	struct Info;

	void setEpsilon(fp_t t) { epsilon = t; }
	void setMaxSubdiv(int v) { maxSubdiv = v; }
};

struct Validator::Subdomain {
	std::vector<Interval> B;
	std::vector<int> qSequence;
	RealInterval time;
	RealInterval incl;

	// Priority function
	bool operator<(const Subdomain &o) const;
	// Debug print
	friend std::ostream& operator<<(std::ostream& ost, const Subdomain &s);
	inline int depth() const { return qSequence.size(); }
	int timeDepth() const {
		int timeDepth = 0;
		const fp_t w = time.width(); 
		if (w == 0) throw std::runtime_error("Zero-width interval");
		fp_t l = 1.;
		while (l > w) {
			l /= 2;
			++timeDepth;
		}
		return timeDepth - depth();
	}
	// Copy the path I took to get here
	void copySequence(std::vector<int> &dst) const;
};


struct Validator::Info {
	int spaceDepth = 0;
	int timeDepth = 0;
	enum class Status {
		unknown, completed, reachedTarget, globalCondition, maxDepth, noSplit
	};
	Status status;

	bool success() const;
	std::string description() const;
};

}