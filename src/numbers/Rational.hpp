#pragma once

#ifdef GMP_INTERFACE

#include <gmpxx.h>

namespace element_validity {
class Rational {
	private:
	mpq_class value;
	
	public:
	// Constructor
	Rational() : value(0) {}
	template <typename T>
	Rational(T v) : value(v) {}
	Rational(int p, int q) : value(p, q) {}

	// Numerator, denominator
	inline mpz_class &num() { return value.get_num(); }
	inline mpz_class &den() { return value.get_den(); }

	// Negation
	inline Rational operator-() const { return -value; }

	// Rat-Rat Comparisons
	inline bool operator>(const Rational &o) const { return value > o.value; }
	inline bool operator>=(const Rational &o) const { return value >= o.value; }
	inline bool operator<(const Rational &o) const { return value < o.value; }
	inline bool operator<=(const Rational &o) const { return value <= o.value; }
	bool operator==(const Rational &o) const { return value == o.value; }
	bool operator!=(const Rational &o) const { return value != o.value; }

	// Rat-FP Comparisons
	inline bool operator>(double o) const { return value > o; }
	inline bool operator>=(double o) const { return value >= o; }
	inline bool operator<(double o) const { return value < o; }
	inline bool operator<=(double o) const { return value <= o; }
	bool operator==(double o) const { return value == o; }
	bool operator!=(double o) const { return value != o; }

	// Rat-Rat Operations
	inline Rational operator+(const Rational &o) const
		{ return value + o.value; }
	inline Rational operator-(const Rational &o) const
		{ return value - o.value; }
	inline Rational operator*(const Rational &o) const
		{ return value * o.value; }
	inline Rational operator/(const Rational &o) const
		{ return value / o.value; }
	inline void operator+=(const Rational &o) { value += o.value; }
	inline void operator-=(const Rational &o) { value -= o.value; }
	inline void operator*=(const Rational &o) { value *= o.value; }
	inline void operator/=(const Rational &o) { value /= o.value; }

	// Rat-FP Operations
	inline Rational operator+(double o) const { return value + o; }
	inline Rational operator-(double o) const { return value - o; }
	inline Rational operator*(double o) const { return value * o; }
	inline Rational operator/(double o) const { return value / o; }
	inline void operator+=(double o) { value += o; }
	inline void operator-=(double o) { value -= o; }
	inline void operator*=(double o) { value *= o; }
	inline void operator/=(double o) { value /= o; }

	// Conversions
	explicit operator double() const { return value.get_d(); }
	explicit operator std::string() const {
		const std::string str = value.get_str();
		// if (str.length() > 8) return std::to_string(den());
		return str;
	}

	friend std::ostream& operator<<(std::ostream& ost, const Rational& r) {
		ost << std::string(r);
		return ost;
	}
};
}

/*
template<> struct Eigen::NumTraits<Rational> :
	Eigen::GenericNumTraits<Rational>
{
	typedef Rational Real;
	typedef Rational NonInteger;
	typedef Rational Nested;
	static inline Real epsilon() { return 0; }
	static inline Real dummy_precision() { return 0; }
	static inline Real digits10() { return 0; }

	enum {
		IsInteger = 0,
		IsSigned = 1,
		IsComplex = 0,
		RequireInitialization = 1,
		ReadCost = 6,
		AddCost = 150,
		MulCost = 100
	};
};
*/

#endif
