#pragma once

#include "../Math/Polynomials.hpp"
#include "../Utility/TypeTraits.hpp"
#include "System.hpp"
#include "../Primitives/Signal.hpp"

#include <vector>


namespace dspbb {


template <class T, eSystemDiscretization Discretization>
class TransferFunctionSystem {
	static_assert(!is_complex_v<T>, "Polynomials must have real coefficients.");

public:
	using C = std::complex<T>;
	TransferFunctionSystem(std::vector<T> numerator = {}, std::vector<T> denominator = {});
	const std::vector<T>& Numerator() const { return m_numerator; }
	const std::vector<T>& Denominator() const { return m_denominator; }
	C operator()(const C& x) const;
	T operator()(const T& x) const;

private:
	std::vector<T> m_numerator;
	std::vector<T> m_denominator;
};

template <class T, eSystemDiscretization Discretization>
TransferFunctionSystem<T, Discretization>::TransferFunctionSystem(std::vector<T> numerator, std::vector<T> denominator)
	: m_numerator(std::move(numerator)), m_denominator(std::move(denominator)) {}


template <class T, eSystemDiscretization Discretization>
typename TransferFunctionSystem<T, Discretization>::C TransferFunctionSystem<T, Discretization>::operator()(const C& x) const {
	return EvaluatePolynomial(m_numerator, x) / EvaluatePolynomial(m_denominator, x);
}

template <class T, eSystemDiscretization Discretization>
T TransferFunctionSystem<T, Discretization>::operator()(const T& x) const {
	return EvaluatePolynomial(m_numerator, x) / EvaluatePolynomial(m_denominator, x);
}


template <class T>
using ContinuousTransferFunctionSystem = TransferFunctionSystem<T, eSystemDiscretization::CONTINUOUS>;

template <class T>
using DiscreteTransferFunctionSystem = TransferFunctionSystem<T, eSystemDiscretization::DISCRETE>;


} // namespace dspbb