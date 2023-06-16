#pragma once

#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Utility/TypeTraits.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <complex>
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
	explicit ApproxComplex(std::complex<double> value) : m_value(value) {}

	template <typename T, typename = typename std::enable_if<std::is_constructible_v<std::complex<double>, T>>::type>
	explicit ApproxComplex(T const& value) : ApproxComplex(static_cast<std::complex<double>>(value)) {}

	template <typename T, typename = typename std::enable_if<std::is_constructible_v<std::complex<double>, T>>::type>
	friend bool operator==(const T& lhs, ApproxComplex const& rhs) {
		auto lhs_v = static_cast<std::complex<double>>(lhs);
		return Catch::Approx(rhs.m_value.real()).epsilon(rhs.m_epsilon).margin(rhs.m_margin).scale(rhs.m_scale) == lhs_v.real()
			   && Catch::Approx(rhs.m_value.imag()).epsilon(rhs.m_epsilon).margin(rhs.m_margin).scale(rhs.m_scale) == lhs_v.imag();
	}

	template <typename T, typename = typename std::enable_if<std::is_constructible_v<std::complex<double>, T>>::type>
	friend bool operator==(ApproxComplex const& lhs, const T& rhs) {
		return operator==(rhs, lhs);
	}

	template <typename T, typename = typename std::enable_if<std::is_constructible_v<std::complex<double>, T>>::type>
	friend bool operator!=(T const& lhs, ApproxComplex const& rhs) {
		return !operator==(lhs, rhs);
	}

	template <typename T, typename = typename std::enable_if<std::is_constructible_v<std::complex<double>, T>>::type>
	friend bool operator!=(ApproxComplex const& lhs, T const& rhs) {
		return !operator==(rhs, lhs);
	}

	template <typename T, typename = typename std::enable_if<std::is_constructible_v<double, T>>::type>
	ApproxComplex& epsilon(T const& newEpsilon) {
		m_epsilon = static_cast<double>(newEpsilon);
		return *this;
	}

	template <typename T, typename = typename std::enable_if<std::is_constructible_v<double, T>>::type>
	ApproxComplex& margin(T const& newMargin) {
		m_margin = static_cast<double>(newMargin);
		return *this;
	}

	template <typename T, typename = typename std::enable_if<std::is_constructible_v<double, T>>::type>
	ApproxComplex& scale(T const& newScale) {
		m_scale = static_cast<double>(newScale);
		return *this;
	}

	std::string toString() const;

	friend std::ostream& operator<<(std::ostream& os, const ApproxComplex& obj) {
		os << obj.m_value;
		return os;
	}

private:
	double m_epsilon = std::numeric_limits<float>::epsilon() * 100.0;
	double m_margin = 0.0;
	double m_scale = 0.0;
	std::complex<double> m_value;
};



template <class T>
dspbb::Signal<T> RandomPositiveSignal(size_t size) {
	thread_local std::mt19937 rne(772537547);
	thread_local std::uniform_real_distribution<float> rng(1, 2);
	dspbb::Signal<T> s;
	for (size_t i = 0; i < size; ++i) {
		s.push_back(rng(rne));
	}
	return s;
}

template <class T, dspbb::eSignalDomain Domain>
dspbb::BasicSignal<T, Domain> RandomSignal(size_t length) {
	thread_local std::mt19937_64 rne(723574);
	thread_local std::uniform_real_distribution<float> rng(-1, 1);
	dspbb::BasicSignal<T, Domain> s(length);
	for (auto& v : s) {
		if constexpr (dspbb::is_complex_v<T>) {
			v.real(rng(rne));
			v.imag(rng(rne));
		}
		else {
			v = rng(rne);
		}
	}
	return s;
}