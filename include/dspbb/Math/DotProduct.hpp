#pragma once

#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Primitives/SignalView.hpp>
#include <dspbb/Vectorization/MathFunctions.hpp>
#include <type_traits>


namespace dspbb {

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU> && !is_complex_v<typename SignalU::value_type>, int> = 0>
auto DotProduct(const SignalT& a, const SignalU& b) {
	assert(a.Size() == b.Size());
	using T = typename SignalT::value_type;
	using U = typename SignalU::value_type;
	using R = decltype(std::declval<T>() * std::declval<U>());
	return InnerProductVectorized(
		a.Data(),
		b.Data(),
		a.Size(),
		R(remove_complex_t<R>(0)),
		[](const auto& a, const auto& b) { return a * b; },
		[](const auto& acc, const auto& x) { return acc + x; });
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU> && is_complex_v<typename SignalU::value_type>, int> = 0>
auto DotProduct(const SignalT& a, const SignalU& b) {
	assert(a.Size() == b.Size());
	using T = typename SignalT::value_type;
	using U = typename SignalU::value_type;
	using R = decltype(std::declval<T>() * std::declval<U>());
	return InnerProductVectorized(
		a.Data(),
		b.Data(),
		a.Size(),
		R(remove_complex_t<R>(0)),
		[](const auto& a, const auto& b) { return a * math_functions::conj(b); },
		[](const auto& acc, const auto& x) { return acc + x; });
}


} // namespace dspbb