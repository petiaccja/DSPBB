#pragma once

#include <catch2/catch.hpp>
#include <complex>
#include <dspbb/Primitives/Signal.hpp>
#include <random>


#define TYPES_COMPLEX          \
	float,                     \
		double,                \
		(std::complex<float>), \
		(std::complex<double>)


#define TYPES_BINARY_REAL ( \
	(float, float),         \
	(float, double),        \
	(double, float),        \
	(double, double))


#define TYPES_BINARY_COMPLEX (                    \
	(float, float),                               \
	(float, double),                              \
	(float, std::complex<float>),                 \
	(double, float),                              \
	(double, double),                             \
	(double, std::complex<double>),               \
	(std::complex<float>, std::complex<float>),   \
	(std::complex<float>, float),                 \
	(std::complex<double>, std::complex<double>), \
	(std::complex<double>, double))

template <class T, class U>
using BinaryProd = decltype(std::declval<T>() * std::declval<U>());
template <class T, class U>
using BinarySum = decltype(std::declval<T>() + std::declval<U>());
template <class T, class U>
using BinaryQuot = decltype(std::declval<T>() / std::declval<U>());
template <class T, class U>
using BinaryDiff = decltype(std::declval<T>() - std::declval<U>());


class ApproxComplex {
public:
	explicit ApproxComplex(std::complex<double> value)
		: m_epsilon(std::numeric_limits<float>::epsilon() * 100),
		  m_margin(0.0),
		  m_scale(0.0),
		  m_value(value) {}

	template <typename T, typename = typename std::enable_if<std::is_constructible<std::complex<double>, T>::value>::type>
	explicit ApproxComplex(T const& value) : ApproxComplex(static_cast<std::complex<double>>(value)) {}

	template <typename T, typename = typename std::enable_if<std::is_constructible<std::complex<double>, T>::value>::type>
	friend bool operator==(const T& lhs, ApproxComplex const& rhs) {
		auto lhs_v = static_cast<std::complex<double>>(lhs);
		return Catch::Detail::Approx(rhs.m_value.real()).epsilon(rhs.m_epsilon).margin(rhs.m_margin).scale(rhs.m_scale) == lhs_v.real()
			   && Catch::Detail::Approx(rhs.m_value.imag()).epsilon(rhs.m_epsilon).margin(rhs.m_margin).scale(rhs.m_scale) == lhs_v.imag();
	}

	template <typename T, typename = typename std::enable_if<std::is_constructible<std::complex<double>, T>::value>::type>
	friend bool operator==(ApproxComplex const& lhs, const T& rhs) {
		return operator==(rhs, lhs);
	}

	template <typename T, typename = typename std::enable_if<std::is_constructible<std::complex<double>, T>::value>::type>
	friend bool operator!=(T const& lhs, ApproxComplex const& rhs) {
		return !operator==(lhs, rhs);
	}

	template <typename T, typename = typename std::enable_if<std::is_constructible<std::complex<double>, T>::value>::type>
	friend bool operator!=(ApproxComplex const& lhs, T const& rhs) {
		return !operator==(rhs, lhs);
	}

	template <typename T, typename = typename std::enable_if<std::is_constructible<double, T>::value>::type>
	ApproxComplex& epsilon(T const& newEpsilon) {
		const double epsilonAsDouble = static_cast<double>(newEpsilon);
		m_epsilon = epsilonAsDouble;
		return *this;
	}

	template <typename T, typename = typename std::enable_if<std::is_constructible<double, T>::value>::type>
	ApproxComplex& margin(T const& newMargin) {
		const double marginAsDouble = static_cast<double>(newMargin);
		m_margin = marginAsDouble;
		return *this;
	}

	template <typename T, typename = typename std::enable_if<std::is_constructible<double, T>::value>::type>
	ApproxComplex& scale(T const& newScale) {
		m_scale = static_cast<double>(newScale);
		return *this;
	}

	std::string toString() const;

private:
	double m_epsilon;
	double m_margin;
	double m_scale;
	std::complex<double> m_value;
};


template <class T>
dspbb::TimeSignal<T> RandomPositiveSignal(size_t size) {
	thread_local std::mt19937 rne(772537547);
	thread_local std::uniform_real_distribution<float> rng(1, 2);
	dspbb::TimeSignal<T> s;
	for (size_t i =0; i<size; ++i) {
		s.PushBack(rng(rne));
	}
	return s;
}