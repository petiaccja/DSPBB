#pragma once

#include "../ComputeKernels/VectorizedAlgorithms.hpp"
#include "../ComputeKernels/VectorizedMathFunctions.hpp"
#include "../Primitives/SignalTraits.hpp"
#include "../Utility/TypeTraits.hpp"

#include <cassert>
#include <type_traits>


namespace dspbb {

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto DotProduct(const SignalT& a, const SignalU& b) {
	assert(a.Size() == b.Size());
	using T = typename SignalT::value_type;
	using U = typename SignalU::value_type;
	using R = product_type_t<T, U>;
	if constexpr (!is_complex_v<U>) {
		return kernels::InnerProductVectorized(
			a.Data(),
			b.Data(),
			a.Size(),
			R(remove_complex_t<R>(0)),
			[](const auto& a, const auto& b) -> product_type_t<decltype(a), decltype(b)> { return a * b; },
			[](const auto& acc, const auto& x) -> sum_type_t<decltype(acc), decltype(x)> { return acc + x; });
	}
	else {
		return kernels::InnerProductVectorized(
			a.Data(),
			b.Data(),
			a.Size(),
			R(remove_complex_t<R>(0)),
			[](const auto& a, const auto& b) -> product_type_t<decltype(a), decltype(b)> { return a * kernels::math_functions::conj(b); },
			[](const auto& acc, const auto& x) -> sum_type_t<decltype(acc), decltype(x)> { return acc + x; });
	}
}

} // namespace dspbb