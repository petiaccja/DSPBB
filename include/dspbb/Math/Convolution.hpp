#pragma once

#include "../Kernels/Convolution.hpp"
#include "../Math/Convolution.hpp"
#include "../Math/DotProduct.hpp"
#include "../Signal/Signal.hpp"
#include "../Signal/SignalView.hpp"
#include "../Utility/TypeTraits.hpp"


namespace dspbb {


enum class eConvolutionMethod {
	/// <summary> Computes only the central part of the convolution, such that the signals don't need to be padded. </summary>
	CENTRAL,
	/// <summary> Computes the full convolution, padding the signals with virtual zeros to both sides. </summary>
	FULL,
};


inline constexpr auto CONV_CENTRAL = std::integral_constant<eConvolutionMethod, eConvolutionMethod::CENTRAL>{};
inline constexpr auto CONV_FULL = std::integral_constant<eConvolutionMethod, eConvolutionMethod::FULL>{};


/// <summary> Calculates the length of the result of the convolution U*V. </summary>
/// <param name="lengthU"> size of U. </param>
/// <param name="lengthV"> size of V. </param>
constexpr size_t ConvolutionLength(size_t lengthU, size_t lengthV, eConvolutionMethod method) {
	if (lengthU == 0 || lengthV == 0) {
		return 0;
	}
	const auto& [shorter, longer] = std::minmax(lengthU, lengthV);
	return method == eConvolutionMethod::CENTRAL ? longer - shorter + 1 : longer + shorter - 1;
}


/// <summary> Calculates the offset to compute central or full convolution. </summary>
/// <param name="lengthU"> size of U. </param>
/// <param name="lengthV"> size of V. </param>
constexpr size_t ConvolutionOffset(size_t lengthU, size_t lengthV, eConvolutionMethod method) {
	if (lengthU == 0 || lengthV == 0) {
		return 0;
	}
	return method == eConvolutionMethod::FULL ? size_t(0) : std::min(lengthU - 1, lengthV - 1);
}


/// <summary> Convolve the two signals. </summary>
/// <param name="out"> The convolved signal is written here. </param>
/// <param name="u"> The first argument to the convolution. </param>
/// <param name="v"> The second argument to the convolution. </param>
/// <param name="offset"> Controls the starting point of the output. </param>
/// <param name="clearOut"> Set to false if the output buffer is already zeroed. </param>
/// <remarks> The subset of the full convolution given by [offset, offset + out.size()).
///		is written into out. </remarks>
template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT, same_domain_as_r<SignalR> SignalU>
void Convolution(SignalR&& out, const SignalT& u, const SignalU& v, size_t offset, bool clearOut = true) {
	const size_t fullLength = ConvolutionLength(u.size(), v.size(), CONV_FULL);
	assert(offset + out.size() <= fullLength && "Result is outside of full convolution, thus contains some true zeros. I mean, it's ok, but you are probably doing it wrong.");

	// Slided is faster, but its accuracy degrades for large input and a compensated reduction is better.
	const size_t shorterSize = std::min(u.size(), v.size());
	if (shorterSize <= 32) {
		kernels::ConvolutionSlide(u.begin(), u.end(), v.begin(), v.end(), out.begin(), out.end(), offset, !clearOut);
	}
	else {
		kernels::ConvolutionReduceVec(u.begin(), u.end(), v.begin(), v.end(), out.begin(), out.end(), offset, !clearOut, plus_compensated<>{});
	}
}


/// <summary> Convolve the two signals. </summary>
/// <param name="out"> The convolved signal is written here. </param>
/// <param name="u"> The first argument to the convolution. </param>
/// <param name="v"> The second argument to the convolution. </param>
/// <param name="clearOut"> Set to false if the output buffer is already zeroed. </param>
template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT, same_domain_as_r<SignalR> SignalU, eConvolutionMethod Method>
void Convolution(SignalR&& out, const SignalT& u, const SignalU& v, std::integral_constant<eConvolutionMethod, Method> method, bool clearOut = true) {
	const auto length = ConvolutionLength(u.size(), v.size(), method);
	const size_t offset = ConvolutionOffset(u.size(), v.size(), method);

	assert(out.size() == length && "Use ConvolutionLength to calculate output size properly.");

	Convolution(out, u, v, offset, clearOut);
}


/// <summary> Convolve the two signals. </summary>
/// <param name="u"> The first argument to the convolution. </param>
/// <param name="v"> The second argument to the convolution. </param>
/// <param name="offset"> Controls the starting point of the output. </param>
/// <param name="length"> Controls the length of the output. </param>
/// <returns> The subset of the full convolution given by [offset, offset + length). </returns>
template <signal_or_view SignalT, same_domain_as<SignalT> SignalU>
auto Convolution(const SignalT& u, const SignalU& v, size_t offset, size_t length) {
	constexpr eSignalDomain Domain = domain_v<std::decay_t<SignalT>>;
	using T = scalar_type_t<std::decay_t<SignalT>>;
	using U = scalar_type_t<std::decay_t<SignalU>>;
	using R = multiplies_result_t<T, U>;

	BasicSignal<R, Domain> out(length, R(remove_complex_t<R>(0)));
	Convolution(out, u, v, offset, false);
	return out;
}


/// <summary> Convolve the two signals. </summary>
/// <param name="u"> The first argument to the convolution. </param>
/// <param name="v"> The second argument to the convolution. </param>
/// <returns> The full or central convolution. </returns>
template <signal_or_view SignalT, same_domain_as<SignalT> SignalU, eConvolutionMethod Method>
auto Convolution(const SignalT& u, const SignalU& v, std::integral_constant<eConvolutionMethod, Method> method) {
	const size_t length = ConvolutionLength(u.size(), v.size(), method);
	const size_t offset = ConvolutionOffset(u.size(), v.size(), method);
	return Convolution(u, v, offset, length);
}

} // namespace dspbb
