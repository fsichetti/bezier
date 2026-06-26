#pragma once

#include <sstream>
#ifdef GMP_INTERFACE
#include "Rational.hpp"
#endif
#include "utils/globals.hpp"
#ifdef IPRED_ARITHMETIC
	#include "numerics.h"
#else
	// #warning Using the naive interval implementation: computations will be slower!
	#include <algorithm>	// for minmax
	#include <cmath>
	#include <stdexcept>
#endif

namespace element_validity::intervals {
	class RobustInterval; // FP Interval with outward rounding operations
	#ifdef EXPERIM_FLOAT_ARITHMETIC
	class FPWrapper; // wrapper to floating point type (NOT ROBUST)
	#endif
}

#ifndef EXPERIM_FLOAT_ARITHMETIC
// Bring names into the element_validity namespace 
namespace element_validity {
	using Interval = intervals::RobustInterval;
	using RealInterval = intervals::RobustInterval;
}
#else
// #warning "You are using non robust arithmetic! Only use this for testing."
namespace element_validity {
	using Interval = intervals::FPWrapper;
	using RealInterval = intervals::RobustInterval;
}
#endif


namespace element_validity::intervals {

#ifdef IPRED_ARITHMETIC
class RobustInterval {
	public:
	static void init();

	private:
	interval_number data;

	inline RobustInterval(const interval_number &o) :
		data(-o.inf(), o.sup()) {}

#ifdef GMP_INTERFACE
	RobustInterval fromRational(const Rational &rat);
#endif

	public:
	// Constructors
	inline RobustInterval(const double &value_lo, const double &value_hi) :
		data(-value_lo, value_hi) {}
	inline RobustInterval(const double &val) : data(val) {}
	inline RobustInterval() : RobustInterval(0.) {}
#ifdef GMP_INTERFACE
	inline RobustInterval(const Rational &rat) : RobustInterval(fromRational(rat)) {}
#endif

	// Access
	inline double lower() const { return data.inf(); }
	inline double upper() const { return data.sup(); }

	// Operators
	inline RobustInterval operator-() const { return -data; }
	inline RobustInterval operator+(const RobustInterval &o) const { return data + o.data; }
	inline RobustInterval operator-(const RobustInterval &o) const { return data - o.data; }
	inline RobustInterval operator*(const RobustInterval &o) const { return data * o.data; }

	inline RobustInterval operator+(double o) const { return data + o; }
	inline friend RobustInterval operator+(double o, const RobustInterval m)
		{ return m + o; }
	inline RobustInterval operator-(double o) const { return data - o; }
	inline friend RobustInterval operator-(double o, const RobustInterval m)
		{ return m - o; }
	inline RobustInterval operator*(double o) const { return data * o; }
	inline friend RobustInterval operator*(double o, const RobustInterval m)
		{ return m * o; }
	inline RobustInterval operator/(double o) const { return data / o; }
	inline friend RobustInterval operator/(double o, const RobustInterval m)
		{ return m / o; }

	inline void operator+=(const RobustInterval &o) { data += o.data; }
	inline void operator-=(const RobustInterval &o) { data -= o.data; }
	inline void operator*=(const RobustInterval &o) { data *= o.data; }

	inline void operator+=(double o) { data += o; }
	inline void operator-=(double o) { data -= o; }
	inline void operator*=(double o) { data *= o; }
	inline void operator/=(double o) { data /= o; }

	inline bool operator<(const RobustInterval &x) const 
		{ return data < x.data; }
	inline bool operator<=(const RobustInterval &x) const 
		{ return data <= x.data; }
	inline bool operator>(const RobustInterval &x) const 
		{ return data > x.data; }
	inline bool operator>=(const RobustInterval &x) const 
		{ return data >= x.data; }
	inline bool operator==(const RobustInterval &x) const 
		{ return data == x.data; }
	inline bool operator!=(const RobustInterval &x) const 
		{ return data != x.data; }

	inline bool operator<(double x) const 
		{ return data < x; }
	inline bool operator<=(double x) const 
		{ return data <= x; }
	inline bool operator>(double x) const 
		{ return data > x; }
	inline bool operator>=(double x) const 
		{ return data >= x; }
	inline bool operator==(double x) const 
		{ return data == x; }
	inline bool operator!=(double x) const 
		{ return data != x; }

	inline RobustInterval abs() const { return data.abs(); }
	inline RobustInterval pow(int e) const { return data.pow(e); }
	inline double width() const { return data.width(); }
	inline double middle() const { return data.getMid(); }
	inline bool contains(double x) const { return (data-x).containsZero(); }
	inline bool contains(const RobustInterval &o) const
		{ return contains(o.lower()) && contains(o.upper()); }

	friend RobustInterval min(const RobustInterval &a, const RobustInterval &b) {
		return min(a.data, b.data);
	}
	friend RobustInterval max(const RobustInterval &a, const RobustInterval &b) {
		return max(a.data, b.data);
	}


	// Conversions
	friend std::ostream& operator<<(std::ostream &ost, const RobustInterval &r);

	explicit operator std::string() const;

	explicit operator double() const { return middle(); }
};

//------------------------------------------------------------------------------

#else

class RobustInterval {
	private:
	fp_t lo, hi;

	public:
	static void init() {};
	// Constructors
	RobustInterval(const fp_t &value_lo, const fp_t &value_hi) :
		lo(value_lo), hi(value_hi) {
		if(lo>hi) throw std::invalid_argument(
			"Lower bound cannot be greater than upper bound");
	}
	RobustInterval(const fp_t &value) : lo(value), hi(value)  {}
	RobustInterval() : lo(0.), hi(0.) {}
	RobustInterval(const RobustInterval &o) : lo(o.lo), hi(o.hi) {}

#ifdef GMP_INTERFACE
	inline RobustInterval(const Rational &rat) : RobustInterval(fromRational(rat)) {}

	RobustInterval fromRational(const Rational &rat);
#endif

	// Access
	inline const fp_t& lower() const { return lo; }
	inline const fp_t& upper() const { return hi; }

	// Opposite
	inline RobustInterval operator-() const { return RobustInterval(-hi, -lo); }

	// Absolute value
	RobustInterval abs() {
		if (hi < 0) return RobustInterval(-hi, -lo);
		if (lo < 0) return RobustInterval(0, std::max(hi, -lo));
		return RobustInterval(lo, hi);
	}

	// Addition
	RobustInterval operator+(const RobustInterval &o) const {
		return RobustInterval(roundDn(lo + o.lo), roundUp(hi + o.hi));
	}
	RobustInterval& operator+=(const RobustInterval& o) {
		lo = roundDn(lo + o.lo);
		hi = roundUp(hi + o.hi);
		return *this;
	}
	RobustInterval operator+(const fp_t &o) const {
		return RobustInterval(roundDn(lo + o), roundUp(hi + o));
	}
	RobustInterval& operator+=(fp_t o) {
		lo = roundDn(lo + o);
		hi = roundUp(hi + o);
		return *this;
	}
	
	// Subtraction
	RobustInterval operator-(const RobustInterval &o) const {
		return RobustInterval(roundDn(lo - o.hi), roundUp(hi - o.lo));
	}
	RobustInterval& operator-=(const RobustInterval& o) {
		lo = roundDn(lo - o.hi);
		hi = roundUp(hi - o.lo);
		return *this;
	}
	RobustInterval operator-(const fp_t &o) const {
		return RobustInterval(roundDn(lo - o), roundUp(hi - o));
	}
	friend RobustInterval operator-(double o, const RobustInterval m) {
		return RobustInterval(roundDn(m.lo - o), roundUp(m.hi - o));
	}
	RobustInterval& operator-=(fp_t o) {
		lo = roundDn(lo - o);
		hi = roundUp(hi - o);
		return *this;
	}

	// Multiplication
	RobustInterval operator*(const RobustInterval &o) const {
		const std::pair<fp_t, fp_t> p = std::minmax({
			lo * o.lo, lo * o.hi,
			hi * o.lo, hi * o.hi
		});
		return RobustInterval(roundDn(p.first), roundUp(p.second));
	}
	RobustInterval& operator*=(const RobustInterval& o) {
		const std::pair<fp_t, fp_t> p = std::minmax({
			lo * o.lo, lo * o.hi,
			hi * o.lo, hi * o.hi
		});
		lo = roundDn(p.first);
		hi = roundUp(p.second);
		return *this;
	}
	RobustInterval operator*(const fp_t &o) const {
		const std::pair<fp_t, fp_t> p = std::minmax({lo * o, hi * o});
		return RobustInterval(roundDn(p.first), roundUp(p.second));
	}
	friend RobustInterval operator*(double o, const RobustInterval m)
	{ 
		const std::pair<fp_t, fp_t> p = std::minmax({m.lo * o, m.hi * o});
		return RobustInterval(roundDn(p.first), roundUp(p.second));
	}
	RobustInterval& operator*=(fp_t o) {
		const std::pair<fp_t, fp_t> p = std::minmax({lo * o, hi * o});
		lo = roundDn(p.first);
		hi = roundUp(p.second);
		return *this;
	}

	// Division by a number
	RobustInterval operator/(fp_t o) const {
		if (o == 0) throw std::runtime_error("Division by 0.");
		const std::pair<fp_t, fp_t> p = std::minmax({lo / o, hi / o});
		return RobustInterval(roundDn(p.first), roundUp(p.second));
	}
	RobustInterval& operator/=(fp_t o) {
		if (o == 0) throw std::runtime_error("Division by 0.");
		const std::pair<fp_t, fp_t> p = std::minmax({lo / o, hi / o});
		lo = roundDn(p.first);
		hi = roundUp(p.second);
		return *this;
	}

	// Square
	RobustInterval sqr() const {
		const std::pair<fp_t, fp_t> p = std::minmax(
			{lo * lo, lo * hi, hi * hi});
		return RobustInterval(roundDn(p.first), roundUp(p.second));
	}

	// Generic power
	RobustInterval pow(int e) const;

	// Whether the interval contains a number or interval
	inline bool contains(fp_t x) const { return (lo <= x && hi >= x); }
	inline bool contains(const RobustInterval &o) const
		{ return contains(o.lower()) && contains(o.upper()); }

	// Comparisons with a floating point number
	inline bool operator<(fp_t x) const { return hi < x; }
	inline bool operator<=(fp_t x) const { return hi <= x; }
	inline bool operator>(fp_t x) const { return lo > x; }
	inline bool operator>=(fp_t x) const { return lo >= x; }
	inline bool operator==(fp_t x) const { return lo == x && hi == x; }
	inline bool operator!=(fp_t x) const { return !contains(x); }

	// Comparisons with another interval
	inline bool operator<(const RobustInterval& o) const { return hi < o.lo; }
	inline bool operator<=(const RobustInterval& o) const { return hi <= o.lo; }
	inline bool operator>(const RobustInterval& o) const { return lo > o.hi; }
	inline bool operator>=(const RobustInterval& o) const { return lo >= o.hi; }
	inline bool operator==(const RobustInterval& o) const
		{ return lo == hi && lo == o.lo && hi == o.hi; }
	inline bool operator!=(const RobustInterval& o) const
		{ return *this < o || *this > o; }

	// Minimum and maximum
	friend inline RobustInterval min(const RobustInterval &a, const RobustInterval &b) {
		return RobustInterval(
			(a.lo < b.lo) ? a.lo : b.lo,
			(a.hi < b.hi) ? a.hi : b.hi
		);
	}
	friend inline RobustInterval max(const RobustInterval &a, const RobustInterval &b) {
		return RobustInterval(
			(a.lo < b.lo) ? b.lo : a.lo,
			(a.hi < b.hi) ? b.hi : a.hi
		);
	}


	// Width or diameter
	inline fp_t width() const { return hi - lo; }

	// Midpoint (approximate)
	inline fp_t middle() const { return lo + (width() / 2); }

	// Conversions
	friend std::ostream& operator<<(std::ostream &ost, const RobustInterval &r);

	explicit operator std::string() const;

	explicit operator fp_t() const { return middle(); }

	private:
	// Rounding up and down
	static constexpr fp_t POSINF = std::numeric_limits<fp_t>::max();
	static constexpr fp_t NEGINF = std::numeric_limits<fp_t>::min();
	inline static fp_t roundUp(const fp_t &value) {
		return std::nextafter(value, POSINF);
	}
	inline static fp_t roundDn(const fp_t &value) {
		return std::nextafter(value, NEGINF);
	}
};
#endif

#ifdef EXPERIM_FLOAT_ARITHMETIC
struct FPWrapper {
	fp_t fp;
	
	static void init() {};
	// Constructors
	FPWrapper() = default;
    FPWrapper(const fp_t &value) : fp(value) {}
    FPWrapper(const FPWrapper &o) = default;
    FPWrapper(FPWrapper&&) = default;
    FPWrapper& operator=(const FPWrapper&) = default;
    FPWrapper& operator=(FPWrapper&&) = default;
    ~FPWrapper() = default;
#ifdef GMP_INTERFACE
	FPWrapper(const Rational &rat) : FPWrapper(fromRational(rat)) {}

	inline FPWrapper fromRational(const Rational &rat)
		{ return static_cast<fp_t>(rat); };
#endif

	// Access
	inline const fp_t& lower() const { return fp; }
	inline const fp_t& upper() const { return fp; }

	// Opposite
	inline FPWrapper operator-() const { return -fp; }

	// Absolute value
	inline FPWrapper abs() { return std::abs(fp); }

	// Addition
	inline FPWrapper operator+(const FPWrapper &o) const { return fp + o.fp; }
	inline FPWrapper& operator+=(const FPWrapper &o) { fp += o.fp; return *this; }
	friend FPWrapper operator+(const fp_t x, const FPWrapper &r) { return r.fp + x; }
	inline FPWrapper& operator+=(fp_t x) { fp += x; return *this; }

	// Subtraction
	inline FPWrapper operator-(const FPWrapper &o) const { return fp - o.fp; }
	inline FPWrapper& operator-=(const FPWrapper &o) { fp -= o.fp; return *this; }
	friend FPWrapper operator-(const fp_t x, const FPWrapper &r) { return r.fp - x; }
	inline FPWrapper& operator-=(fp_t x) { fp -= x; return *this; }

	// Multiplication
	inline FPWrapper operator*(const FPWrapper &o) const { return fp * o.fp; }
	inline FPWrapper& operator*=(const FPWrapper &o) { fp *= o.fp; return *this; }
	friend FPWrapper operator*(const fp_t x, const FPWrapper &r) { return r.fp * x; }
	inline FPWrapper& operator*=(fp_t x) { fp *= x; return *this; }


	inline FPWrapper operator/(fp_t x) { return fp / x; }
	inline FPWrapper& operator/=(fp_t x) { fp /= x; return *this; }

	inline FPWrapper sqr(int e) const { return fp*fp; };
	inline FPWrapper pow(int e) const { return std::pow(fp, e); };

	// Whether the interval contains a number or interval
	// inline bool contains(fp_t x) const { return (lo <= x && hi >= x); }
	// inline bool contains(const FPWrapper &o) const { return fp == o.fp; }

	// Comparisons with a floating point number
	inline bool operator<(fp_t x) const { return fp < x; }
	inline bool operator<=(fp_t x) const { return fp <= x; }
	inline bool operator>(fp_t x) const { return fp > x; }
	inline bool operator>=(fp_t x) const { return fp >= x; }
	inline bool operator==(fp_t x) const { return fp == x; }
	inline bool operator!=(fp_t x) const { return fp != x; }

	// Comparisons with another interval
	inline bool operator<(const FPWrapper& o) const { return fp < o.fp; }
	inline bool operator<=(const FPWrapper& o) const { return fp <= o.fp; }
	inline bool operator>(const FPWrapper& o) const { return fp > o.fp; }
	inline bool operator>=(const FPWrapper& o) const { return fp >= o.fp; }
	inline bool operator==(const FPWrapper& o) const { return fp == o.fp; }
	inline bool operator!=(const FPWrapper& o) const { return fp != o.fp; }

	// Minimum and maximum
	friend inline FPWrapper min(const FPWrapper &a, const FPWrapper &b) {
		return std::min(a.fp, b.fp);
	}
	friend inline FPWrapper max(const FPWrapper &a, const FPWrapper &b) {
		return std::max(a.fp, b.fp);
	}


	// Width or diameter
	// inline fp_t width() const { return hi - lo; }

	// Midpoint (approximate)
	inline fp_t middle() const { return fp; }

	// Conversions
	friend std::ostream& operator<<(std::ostream &ost, const FPWrapper &r) {
		return ost << r;
	}

	explicit operator std::string() const { return std::to_string(fp); }

	explicit operator fp_t() const { return fp; }
};
#endif

}

/*
template<> struct Eigen::NumTraits<RobustInterval> :
	Eigen::GenericNumTraits<RobustInterval>
{
	typedef RobustInterval Real;
	typedef RobustInterval NonInteger;
	typedef RobustInterval Nested;
	static inline Real epsilon() { return 0; }
	static inline Real dummy_precision() { return 0; }
	static inline Real digits10() { return 0; }

	enum {
		IsInteger = 0,
		IsSigned = 1,
		IsComplex = 0,
		RequireInitialization = 1,
		ReadCost = 2,
		AddCost = 4,
		MulCost = 30
	};
};
*/