#pragma once

#include "../Math/Polynomials.hpp"


namespace dspbb {

enum class eSystemDiscretization {
	DISCRETE,
	CONTINUOUS,
};


template <class T, eSystemDiscretization Discretization>
struct PoleZeroSystem {
	T gain;
	FactoredPolynomial<T> zeros;
	FactoredPolynomial<T> poles;

	std::complex<T> operator()(const std::complex<T>& x) const { return gain * zeros(x) / poles(x); }
	T operator()(const T& x) const { return gain * zeros(x) / poles(x); }
};


template <class T, eSystemDiscretization Discretization>
struct TransferFunctionSystem {
	Polynomial<T> numerator;
	Polynomial<T> denominator;

	std::complex<T> operator()(const std::complex<T>& x) const { return numerator(x) / denominator(x); }
	T operator()(const T& x) const { return numerator(x) / denominator(x); }
};


template <class T>
using ContinuousTransferFunctionSystem = TransferFunctionSystem<T, eSystemDiscretization::CONTINUOUS>;
template <class T>
using DiscreteTransferFunctionSystem = TransferFunctionSystem<T, eSystemDiscretization::DISCRETE>;

template <class T>
using ContinuousPoleZeroSystem = PoleZeroSystem<T, eSystemDiscretization::CONTINUOUS>;
template <class T>
using DiscretePoleZeroSystem = PoleZeroSystem<T, eSystemDiscretization::DISCRETE>;


} // namespace dspbb