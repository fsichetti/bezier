#include "Interval.hpp"

namespace element_validity::intervals {

#ifdef IPRED_ARITHMETIC

void RobustInterval::init() {
	#ifdef IPRED_ARITHMETIC
	static bool initialized = false;
	if (initialized) return;
	initFPU();
	setFPUModeToRoundUP();
	initialized = true;
	#endif
};

#ifdef GMP_INTERFACE
RobustInterval RobustInterval::fromRational(const Rational &rat) {
	// TODO change this function to not use nextafter
	const double inf = std::numeric_limits<double>::max();
	const double d = static_cast<double>(rat);

	if (rat < 0) return RobustInterval(std::nextafter(d, -inf), d);
	if (rat > 0) return RobustInterval(d, std::nextafter(d, inf));
	return RobustInterval(0);
}
#endif

#else

RobustInterval RobustInterval::pow(int e) const {
	switch (e) {
	// Base cases
	case 0: return RobustInterval(1);
	case 1: return *this;
	case 2: return sqr();
	default:
		// Exponent is even: divide by 2 and recurse
		if ((e % 2) == 0) return pow(e >> 1).sqr();
		// Exponent is not divisible by 2
		return *this * pow(e - 1);
	}
}

#ifdef GMP_INTERFACE
RobustInterval RobustInterval::fromRational(const Rational &rat) {
		const double d = static_cast<double>(rat);
		if (rat < 0) return RobustInterval(roundDn(d), d);
		if (rat > 0) return RobustInterval(d, roundUp(d));
		return RobustInterval(0);
	}
#endif
#endif

std::ostream& operator<<(std::ostream &ost, const RobustInterval &r) {
	const auto p = std::setprecision(
		std::numeric_limits<double>::max_digits10);
	ost << '[' << p << r.lower() << ',' << p << r.upper() << ']';
	return ost;
}

RobustInterval::operator std::string() const {
	std::stringstream s;
	s << *this;
	return s.str();
}

}
