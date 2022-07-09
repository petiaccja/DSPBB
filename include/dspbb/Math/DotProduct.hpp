#pragma once

#include "../Kernels/Math.hpp"
#include "../Kernels/Numeric.hpp"
#include "../Primitives/SignalTraits.hpp"
#include "../Utility/TypeTraits.hpp"

#include <cassert>
#include <type_traits>


namespace dspbb {

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto DotProduct(const SignalT& a, const SignalU& b) {
	assert(a.size() == b.size());
	using T = typename SignalT::value_type;
	using U = typename SignalU::value_type;
	using R = multiplies_result_t<T, U>;
	if constexpr (!is_complex_v<U>) {
		return kernels::InnerProduct(
			a.begin(),
			a.end(),
			b.begin(),
			R(remove_complex_t<R>(0)),
			[](const auto& acc, const auto& x) -> plus_result_t<decltype(acc), decltype(x)> { return acc + x; },
			[](const auto& a, const auto& b) -> multiplies_result_t<decltype(a), decltype(b)> { return a * b; });
	}
	else {
		return kernels::InnerProduct(
			a.begin(),
			a.end(),
			b.begin(),
			R(remove_complex_t<R>(0)),
			[](const auto& acc, const auto& x) -> plus_result_t<decltype(acc), decltype(x)> { return acc + x; },
			[](const auto& a, const auto& b) -> multiplies_result_t<decltype(a), decltype(b)> { return a * kernels::math_functions::conj(b); });
	}
}

} // namespace dspbb