#include "Validator.hpp"

namespace element_validity {

RealInterval Validator::minclusion(
	const std::vector<Interval>& B
) const {
	Interval lo(std::numeric_limits<fp_t>::max());
	for (const Interval &b : B) lo = min(lo, b);
	Interval hi(std::numeric_limits<fp_t>::max());
	for (const int c : interpIndices) hi = min(hi, B.at(c));
	return {lo.lower(), hi.upper()};
}

// Priority function
bool Validator::Subdomain::operator<(const Validator::Subdomain &o) const {
	if (time.lower() != o.time.lower())
		return time.lower() > o.time.lower();
	else return incl.lower() > o.incl.lower();
}


// Debug print
std::ostream& operator<<(std::ostream& ost, const Validator::Subdomain &s) {
	ost << "t: " << s.time << std::endl;
	ost << "I: " << s.incl << std::endl;
	ost << "Q: ";
	for (const int x : s.qSequence) ost << x << ", ";
	ost << std::endl;
	ost << "B: ";
	for (const Interval &b : s.B) ost << '\t' << b << std::endl;

	return ost;
}

// Copy the path I took to get here
void Validator::Subdomain::copySequence(std::vector<int> &dst) const {
	dst.clear();
	dst.reserve(depth());
	std::copy(
		qSequence.cbegin(), qSequence.cend(), std::back_inserter(dst)
	);
}

bool Validator::Info::success() const {
	switch (status) {
	case Status::completed:
	case Status::reachedTarget:
	case Status::globalCondition:
		return true;
	default:
		return false;
	}
}
std::string Validator::Info::description() const {
	switch (status) {
	case Status::completed:
		return "Processed all intervals";
	case Status::reachedTarget:
		return "Reached target precision";
	case Status::maxDepth:
		return "User specified max depth reached";
	case Status::globalCondition:
		return "Global early termination condition satisfied";
	case Status::noSplit:
		return "Cannot split due to machine precision";
	default: return "Something is wrong";
	}
}

}