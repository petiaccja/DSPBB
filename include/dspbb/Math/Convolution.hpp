#pragma once

#include "../Kernels/Convolution.hpp"
#include "../Math/Convolution.hpp"
#include "../Math/DotProduct.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "../Utility/TypeTraits.hpp"

#include <complex>

namespace dspbb {

namespace impl {
	class ConvCentral {};
	class ConvFull {};
	constexpr ConvCentral CONV_CENTRAL;
	constexpr ConvFull CONV_FULL;
} // namespace impl

using impl::CONV_CENTRAL;
using impl::CONV_FULL;

/// <summary> Calculates the length of the result of the convolution U*V. </summary>
/// <param name="lengthU"> Length of U. </param>
/// <param name="lengthV"> Length of V. </param>
inline size_t ConvolutionLength(size_t lengthU, size_t lengthV, impl::ConvCentral) {
	if (lengthU == 0 || lengthV == 0) {
		return 0;
	}
	const auto& mm = std::minmax(lengthU, lengthV);
	const auto& shorter = mm.first;
	const auto& longer = mm.second;
	return longer - shorter + 1;
}

/// <summary> Calculates the length of the result of the convolution U*V. </summary>
/// <param name="lengthU"> Length of U. </param>
/// <param name="lengthV"> Length of V. </param>
inline size_t ConvolutionLength(size_t lengthU, size_t lengthV, impl::ConvFull) {
	if (lengthU == 0 || lengthV == 0) {
		return 0;
	}
	const auto& mm = std::minmax(lengthU, lengthV);
	const auto& shorter = mm.first;
	const auto& longer = mm.second;
	return longer + shorter - 1;
}

template <class SignalR, class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalR, SignalT, SignalU>, int> = 0>
auto Convolution(SignalR&& out, const SignalT& u, const SignalU& v, size_t offset, bool clearOut = true) {
	const size_t fullLength = ConvolutionLength(u.Length(), v.Length(), CONV_FULL);
	assert(offset + out.Size() <= fullLength && "Result is outside of full convolution, thus contains some true zeros. I mean, it's ok, but you are probably doing it wrong.");

	// Slided is faster, but it's accuracy degrades for large input and a compensated reduction is better.
	const size_t shorterSize = std::min(u.Length(), v.Length());
	if (shorterSize <= 32) {
		kernels::ConvolutionSlide(u.begin(), u.end(), v.begin(), v.end(), out.begin(), out.end(), offset, !clearOut);
	}
	else {
		kernels::ConvolutionReduceVec(u.begin(), u.end(), v.begin(), v.end(), out.begin(), out.end(), offset, !clearOut, plus_compensated<>{});
	}
}

template <class SignalR, class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalR, SignalT, SignalU>, int> = 0>
auto Convolution(SignalR&& out, const SignalT& u, const SignalU& v, impl::ConvFull, bool clearOut = true) {
	const size_t fullLength = ConvolutionLength(u.Length(), v.Length(), CONV_FULL);
	assert(out.Size() == fullLength && "Use ConvolutionLength to calculate output size properly.");
	const size_t offset = 0;
	Convolution(out, u, v, offset, clearOut);
}

template <class SignalR, class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalR, SignalT, SignalU>, int> = 0>
auto Convolution(SignalR&& out, const SignalT& u, const SignalU& v, impl::ConvCentral, bool clearOut = true) {
	const size_t centralLength = ConvolutionLength(u.Length(), v.Length(), CONV_CENTRAL);
	assert(out.Size() == centralLength && "Use ConvolutionLength to calculate output size properly.");
	const size_t offset = std::min(u.Size() - 1, v.Size() - 1);
	Convolution(out, u, v, offset, clearOut);
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto Convolution(const SignalT& u, const SignalU& v, size_t offset, size_t length) {
	constexpr eSignalDomain Domain = signal_traits<std::decay_t<SignalT>>::domain;
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	using U = typename signal_traits<std::decay_t<SignalU>>::type;
	using R = multiplies_result_t<T, U>;

	BasicSignal<R, Domain> out(length, R(0));
	Convolution(out, u, v, offset, false);
	return out;
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto Convolution(const SignalT& u, const SignalU& v, impl::ConvFull) {
	const size_t length = ConvolutionLength(u.Length(), v.Length(), CONV_FULL);
	const size_t offset = 0;
	return Convolution(u, v, offset, length);
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto Convolution(const SignalT& u, const SignalU& v, impl::ConvCentral) {
	const size_t length = ConvolutionLength(u.Length(), v.Length(), CONV_CENTRAL);
	const size_t offset = std::min(u.Size() - 1, v.Size() - 1);
	return Convolution(u, v, offset, length);
}

} // namespace dspbb
