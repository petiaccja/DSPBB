#pragma once

#include "../Math/Polynomials.hpp"


namespace dspbb {

enum class eSystemDiscretization {
	DISCRETE,
	CONTINUOUS,
};


template <class T, eSystemDiscretization Discretization>
struct ZeroPoleGain {
	T gain;
	FactoredPolynomial<T> zeros;
	FactoredPolynomial<T> poles;

	std::complex<T> operator()(const std::complex<T>& x) const { return gain * zeros(x) / poles(x); }
	T operator()(const T& x) const { return gain * zeros(x) / poles(x); }
};


template <class T, eSystemDiscretization Discretization>
struct TransferFunction {
	TransferFunction() = default;
	TransferFunction(const ZeroPoleGain<T, Discretization>& zpk);

	Polynomial<T> numerator;
	Polynomial<T> denominator;

	std::complex<T> operator()(const std::complex<T>& x) const { return numerator(x) / denominator(x); }
	T operator()(const T& x) const { return numerator(x) / denominator(x); }
};


template <class T>
struct CascadedBiquad {
	struct Biquad {
		std::array<T, 3> numerator;
		std::array<T, 2> denominator;
	};
	std::vector<Biquad> sections;
};


template <class T, eSystemDiscretization Discretization>
TransferFunction<T, Discretization>::TransferFunction(const ZeroPoleGain<T, Discretization>& zpk)
	: numerator{ ExpandPolynomial(zpk.zeros) },
	  denominator{ ExpandPolynomial(zpk.poles) } {
	numerator.Coefficients() *= zpk.gain;
}



template <class T>
using ContinuousTransferFunction = TransferFunction<T, eSystemDiscretization::CONTINUOUS>;
template <class T>
using DiscreteTransferFunction = TransferFunction<T, eSystemDiscretization::DISCRETE>;

template <class T>
using ContinuousZeroPoleGain = ZeroPoleGain<T, eSystemDiscretization::CONTINUOUS>;
template <class T>
using DiscreteZeroPoleGain = ZeroPoleGain<T, eSystemDiscretization::DISCRETE>;


} // namespace dspbb