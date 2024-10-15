#pragma once

#include "../../Math/Convolution.hpp"
#include "../../Math/OverlapAdd.hpp"
#include "../../Signal/Traits.hpp"
#include "../../Utility/TypeTraits.hpp"

#include <cassert>

namespace dspbb {


enum class eFilterMethod {
	/// <summary> Apply FIR filter using convolution. </summary>
	CONVOLUTION,

	/// <summary> Apply FIR filter using the ovarlap-add method via FFTs. </summary>
	OVERLAP_ADD,
};


inline constexpr auto FILTER_CONV = std::integral_constant<eFilterMethod, eFilterMethod::CONVOLUTION>{};
inline constexpr auto FILTER_OLA = std::integral_constant<eFilterMethod, eFilterMethod::OVERLAP_ADD>{};


namespace impl {

	template <mutable_signal_or_view_r SignalS, same_domain_as_r<SignalS> SignalU>
	void ShiftFilterState(SignalS&& state, const SignalU& signal) {
		if (signal.size() < state.size()) {
			std::move(state.begin() + signal.size(), state.end(), state.begin());
		}
		std::copy(signal.rbegin(), signal.rbegin() + std::min(signal.size(), state.size()), state.rbegin());
	}

	template <signal_or_view SignalT, same_domain_as<SignalT> SignalU>
	using ProductSignalT = BasicSignal<multiplies_result_t<typename std::decay_t<SignalT>::value_type, typename std::decay_t<SignalU>::value_type>, domain_v<std::decay_t<SignalT>>>;

} // namespace impl



/// <summary> Calculate the size of the state required for chunk-based FFT filtering. </summary>
/// <param name="filterSize"> The size of the FIR filter. </param>
constexpr size_t FilterStateSize(size_t filterSize) {
	return filterSize - 1;
}


/// <summary> Apply an FIR filter to a signal. </summary>
/// <param name="out"> Output buffer for the filtered signal. </param>
/// <param name="signal"> The signal to filter. </param>
/// <param name="filter"> The filter to apply. </param>
/// <param name="chunkSize"> The FFT's size if the overlap-add method is used. </param>
template <mutable_signal_or_view_r SignalR,
		  same_domain_as_r<SignalR> SignalU,
		  same_domain_as_r<SignalR> SignalV,
		  eConvolutionMethod ConvMethod,
		  eFilterMethod FilterMethod>
void Filter(SignalR&& out,
			const SignalU& signal,
			const SignalV& filter,
			std::integral_constant<eConvolutionMethod, ConvMethod> convMethod,
			std::integral_constant<eFilterMethod, FilterMethod>,
			[[maybe_unused]] size_t chunkSize = 0) {
	if constexpr (FilterMethod == eFilterMethod::CONVOLUTION) {
		Convolution(out, signal, filter, convMethod);
	}
	else {
		OverlapAdd(out, signal, filter, convMethod, chunkSize);
	}
}


/// <summary> Apply an FIR filter to a chunk of a signal, saving state for the next chunk. </summary>
/// <param name="out"> Output buffer for the filtered signal. </param>
/// <param name="signal"> The signal to filter. </param>
/// <param name="filter"> The filter to apply. </param>
/// <param name="state"> State carried over to the next filter call when processing in batches. </param>
/// <param name="chunkSize"> The FFT's size if the overlap-add method is used. </param>
template <mutable_signal_or_view_r SignalR,
		  same_domain_as_r<SignalR> SignalU,
		  same_domain_as_r<SignalR> SignalV,
		  mutable_signal_or_view_r SignalS,
		  eFilterMethod FilterMethod>
void Filter(SignalR&& out,
			const SignalU& signal,
			const SignalV& filter,
			SignalS&& state,
			std::integral_constant<eFilterMethod, FilterMethod>,
			size_t chunkSize = 0) {
	assert(state.size() == FilterStateSize(filter.size()));
	assert(out.size() == signal.size());

	std::fill(out.begin(), out.end(), remove_complex_t<typename std::decay_t<SignalR>::value_type>(0));
	const auto outHead = AsView(out).subsignal(0, std::min(out.size(), state.size()));
	if constexpr (FilterMethod == eFilterMethod::CONVOLUTION) {
		Convolution(outHead, state, filter, filter.size() - 1, false);
		Convolution(out, signal, filter, 0, false);
	}
	else {
		OverlapAdd(outHead, state, filter, filter.size() - 1, chunkSize, false);
		OverlapAdd(out, signal, filter, 0, chunkSize, false);
	}
	impl::ShiftFilterState(state, signal);
}


/// <summary> Apply an FIR filter to a signal. </summary>
/// <param name="signal"> The signal to filter. </param>
/// <param name="filter"> The filter to apply. </param>
/// <param name="chunkSize"> The FFT's size if the overlap-add method is used. </param>
template <signal_or_view SignalU,
		  same_domain_as<SignalU> SignalV,
		  eConvolutionMethod ConvMethod,
		  eFilterMethod FilterMethod>
auto Filter(const SignalU& signal,
			const SignalV& filter,
			std::integral_constant<eConvolutionMethod, ConvMethod> convMethod,
			std::integral_constant<eFilterMethod, FilterMethod> filterMethod,
			size_t chunkSize = 0) {
	impl::ProductSignalT<SignalU, SignalV> out(ConvolutionLength(signal.size(), filter.size(), ConvMethod));
	Filter(out, signal, filter, convMethod, filterMethod, chunkSize);
	return out;
}


/// <summary> Apply an FIR filter to a chunk of a signal, saving state for the next chunk. </summary>
/// <param name="signal"> The signal to filter. </param>
/// <param name="filter"> The filter to apply. </param>
/// <param name="state"> State carried over to the next filter call when processing in batches. </param>
/// <param name="chunkSize"> The FFT's size if the overlap-add method is used. </param>
template <signal_or_view SignalU,
		  same_domain_as<SignalU> SignalV,
		  mutable_signal_or_view_r SignalS,
		  eFilterMethod FilterMethod>
auto Filter(const SignalU& signal,
			const SignalV& filter,
			SignalS&& state,
			std::integral_constant<eFilterMethod, FilterMethod> filterMethod,
			size_t chunkSize = 0) {
	impl::ProductSignalT<SignalU, SignalV> out(signal.size());
	Filter(out, signal, filter, state, filterMethod, chunkSize);
	return out;
}

} // namespace dspbb