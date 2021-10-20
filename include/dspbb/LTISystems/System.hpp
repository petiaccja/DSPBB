#pragma once

#include "../Math/Polynomials.hpp"


namespace dspbb {

enum class eSystemDiscretization {
	DISCRETE,
	CONTINUOUS,
};


template <class T, eSystemDiscretization Discretization>
struct PoleZeroSystem_2 {
	T gain;
	FactoredPolynomial<T> zeros;
	FactoredPolynomial<T> poles;

	std::complex<T> operator()(const std::complex<T>& x) const { return gain * zeros(x) / poles(x); }
	T operator()(const T& x) const { return gain * zeros(x) / poles(x); }
};


template <class T, eSystemDiscretization Discretization>
struct TransferFunctionSystem_2 {
	Polynomial<T> numerator;
	Polynomial<T> denominator;

	std::complex<T> operator()(const std::complex<T>& x) const { return numerator(x) / denominator(x); }
	T operator()(const T& x) const { return numerator(x) / denominator(x); }
};


} // namespace dspbb