#pragma once

#include "../Math/Polynomials.hpp"
#include "System.hpp"

#include <algorithm>


namespace dspbb {

template <class T, eSystemDiscretization Discretization>
TransferFunctionSystem<T, Discretization> TransferFunction(const PoleZeroSystem<T, Discretization>& pz) {
	Polynomial<T> num = ExpandPolynomial(pz.zeros);
	Polynomial<T> den = ExpandPolynomial(pz.poles);
	num.Coefficients() *= pz.gain;
	return { std::move(num), std::move(den) };
}

} // namespace dspbb