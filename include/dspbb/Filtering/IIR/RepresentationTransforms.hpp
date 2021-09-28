#pragma once

#include "../../Math/Polynomials.hpp"
#include "PoleZeroSystem.hpp"
#include "TransferFunctionSystem.hpp"

#include <algorithm>


namespace dspbb {

template <class T, eSystemDiscretization Discretization>
TransferFunctionSystem<T, Discretization> TransferFunction(const PoleZeroSystem<T, Discretization>& pz) {
	std::vector<T> num = ExpandPolynomialReal(pz.Zeros());
	std::vector<T> den = ExpandPolynomialReal(pz.Poles());
	std::for_each(num.begin(), num.end(), [gain = pz.Gain()](auto& arg) { arg *= gain; });
	return { std::move(num), std::move(den) };
}

} // namespace dspbb