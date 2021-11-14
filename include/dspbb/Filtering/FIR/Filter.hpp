#pragma once

#include "../../Math/Convolution.hpp"
#include "../../Math/OverlapAdd.hpp"
#include "../../Primitives/SignalTraits.hpp"
#include "../../Utility/TypeTraits.hpp"


namespace dspbb {

namespace impl {
	struct FilterConv {};
	struct FilterOla {};
	constexpr FilterConv FILTER_CONV;
	constexpr FilterOla FILTER_OLA;


	template <class SignalS, class SignalU>
	void ShiftFilterState(SignalS& state, const SignalU& signal) {
		if (signal.Size() < state.Size()) {
			std::move(state.begin() + signal.Size(), state.end(), state.begin());
		}
		std::copy(signal.rbegin(), signal.rbegin() + std::min(signal.Size(), state.Size()), state.rbegin());
	}

	template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
	using ProductSignalT = Signal<product_type_t<typename std::decay_t<SignalT>::value_type, typename std::decay_t<SignalU>::value_type>, signal_traits<std::decay_t<SignalT>>::domain>;
} // namespace impl


using impl::FILTER_CONV;
using impl::FILTER_OLA;


template <class SignalR, class SignalU, class SignalV, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalU, SignalV>, int> = 0>
auto Filter(SignalR&& out, const SignalU& signal, const SignalV& filter, impl::ConvCentral, impl::FilterOla, size_t chunkSize) {
	OverlapAdd(out, signal, filter, CONV_CENTRAL, chunkSize);
}

template <class SignalR, class SignalU, class SignalV, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalU, SignalV>, int> = 0>
auto Filter(SignalR&& out, const SignalU& signal, const SignalV& filter, impl::ConvCentral, impl::FilterConv) {
	Convolution(out, signal, filter, CONV_CENTRAL);
}

template <class SignalR, class SignalU, class SignalV, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalU, SignalV>, int> = 0>
auto Filter(SignalR&& out, const SignalU& signal, const SignalV& filter, impl::ConvFull, impl::FilterOla, size_t chunkSize) {
	OverlapAdd(out, signal, filter, CONV_FULL, chunkSize);
}

template <class SignalR, class SignalU, class SignalV, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalU, SignalV>, int> = 0>
auto Filter(SignalR&& out, const SignalU& signal, const SignalV& filter, impl::ConvFull, impl::FilterConv) {
	Convolution(out, signal, filter, CONV_FULL);
}

template <class SignalR,
		  class SignalU,
		  class SignalV,
		  class SignalS,
		  std::enable_if_t<is_mutable_signal_v<SignalR> && is_mutable_signal_v<SignalS> && is_same_domain_v<SignalR, SignalU, SignalV, SignalS>, int> = 0>
auto Filter(SignalR&& out, const SignalU& signal, const SignalV& filter, SignalS& state, impl::FilterOla, size_t chunkSize) {
	if (state.Size() != filter.Size() - 1) {
		throw std::invalid_argument("State must have a length one less than the filter.");
	}
	if (out.Size() != signal.Size()) {
		throw std::invalid_argument("Output and input signals must have the same size.");
	}

	std::fill(out.begin(), out.end(), remove_complex_t<typename std::decay_t<SignalR>::value_type>(0));
	OverlapAdd(AsView(out).SubSignal(0, std::min(out.Size(), state.Size())), state, filter, filter.Size() - 1, chunkSize, false);
	OverlapAdd(out, signal, filter, 0, chunkSize, false);
	impl::ShiftFilterState(state, signal);
}

template <class SignalR,
		  class SignalU,
		  class SignalV,
		  class SignalS,
		  std::enable_if_t<is_mutable_signal_v<SignalR> && is_mutable_signal_v<SignalS> && is_same_domain_v<SignalR, SignalU, SignalV, SignalS>, int> = 0>
auto Filter(SignalR&& out, const SignalU& signal, const SignalV& filter, SignalS& state, impl::FilterConv) {
	if (state.Size() != filter.Size() - 1) {
		throw std::invalid_argument("State must have a length one less than the filter.");
	}
	if (out.Size() != signal.Size()) {
		throw std::invalid_argument("Output and input signals must have the same size.");
	}

	std::fill(out.begin(), out.end(), remove_complex_t<typename std::decay_t<SignalR>::value_type>(0));
	Convolution(AsView(out).SubSignal(0, std::min(out.Size(), state.Size())), state, filter, filter.Size() - 1, false);
	Convolution(out, signal, filter, 0, false);
	impl::ShiftFilterState(state, signal);
}

template <class SignalU, class SignalV, std::enable_if_t<is_same_domain_v<SignalU, SignalV>, int> = 0>
auto Filter(const SignalU& signal, const SignalV& filter, impl::ConvCentral, impl::FilterOla, size_t chunkSize) {
	impl::ProductSignalT<SignalU, SignalV> out(ConvolutionLength(signal.Size(), filter.Size(), CONV_CENTRAL));
	Filter(out, signal, filter, CONV_CENTRAL, FILTER_OLA, chunkSize);
	return out;
}

template <class SignalU, class SignalV, std::enable_if_t<is_same_domain_v<SignalU, SignalV>, int> = 0>
auto Filter(const SignalU& signal, const SignalV& filter, impl::ConvCentral, impl::FilterConv) {
	impl::ProductSignalT<SignalU, SignalV> out(ConvolutionLength(signal.Size(), filter.Size(), CONV_CENTRAL));
	Filter(out, signal, filter, CONV_CENTRAL, FILTER_CONV);
	return out;
}

template <class SignalU, class SignalV, std::enable_if_t<is_same_domain_v<SignalU, SignalV>, int> = 0>
auto Filter(const SignalU& signal, const SignalV& filter, impl::ConvFull, impl::FilterOla, size_t chunkSize) {
	impl::ProductSignalT<SignalU, SignalV> out(ConvolutionLength(signal.Size(), filter.Size(), CONV_FULL));
	Filter(out, signal, filter, CONV_FULL, FILTER_OLA, chunkSize);
	return out;
}

template <class SignalU, class SignalV, std::enable_if_t<is_same_domain_v<SignalU, SignalV>, int> = 0>
auto Filter(const SignalU& signal, const SignalV& filter, impl::ConvFull, impl::FilterConv) {
	impl::ProductSignalT<SignalU, SignalV> out(ConvolutionLength(signal.Size(), filter.Size(), CONV_FULL));
	Filter(out, signal, filter, CONV_FULL, FILTER_CONV);
	return out;
}

template <class SignalU,
		  class SignalV,
		  class SignalS,
		  std::enable_if_t<is_mutable_signal_v<SignalS> && is_same_domain_v<SignalU, SignalV, SignalS>, int> = 0>
auto Filter(const SignalU& signal, const SignalV& filter, SignalS&& state, impl::FilterOla, size_t chunkSize) {
	impl::ProductSignalT<SignalU, SignalV> out(signal.Size());
	Filter(out, signal, filter, state, FILTER_OLA, chunkSize);
	return out;
}

template <class SignalU,
		  class SignalV,
		  class SignalS,
		  std::enable_if_t<is_mutable_signal_v<SignalS> && is_same_domain_v<SignalU, SignalV, SignalS>, int> = 0>
auto Filter(const SignalU& signal, const SignalV& filter, SignalS&& state, impl::FilterConv) {
	impl::ProductSignalT<SignalU, SignalV> out(signal.Size());
	Filter(out, signal, filter, state, FILTER_CONV);
	return out;
}

} // namespace dspbb