#pragma once

#include "../../Math/Convolution.hpp"
#include "../../Math/OverlapAdd.hpp"
#include "../../Primitives/SignalTraits.hpp"
#include "../../Utility/TypeTraits.hpp"

#include <cassert>

namespace dspbb {

namespace impl {
	struct FilterConv {};
	struct FilterOla {};
	constexpr FilterConv FILTER_CONV;
	constexpr FilterOla FILTER_OLA;


	template <class SignalS, class SignalU>
	void ShiftFilterState(SignalS& state, const SignalU& signal) {
		if (signal.size() < state.size()) {
			std::move(state.begin() + signal.size(), state.end(), state.begin());
		}
		std::copy(signal.rbegin(), signal.rbegin() + std::min(signal.size(), state.size()), state.rbegin());
	}

	template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
	using ProductSignalT = BasicSignal<multiplies_result_t<typename std::decay_t<SignalT>::value_type, typename std::decay_t<SignalU>::value_type>, signal_traits<std::decay_t<SignalT>>::domain>;
} // namespace impl


using impl::FILTER_CONV;
using impl::FILTER_OLA;


template <class SignalR, class SignalU, class SignalV, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalU, SignalV>, int> = 0>
auto Filter(SignalR&& out, const SignalU& signal, const SignalV& filter, impl::ConvCentral, impl::FilterOla, size_t chunkSize = 0) {
	OverlapAdd(out, signal, filter, CONV_CENTRAL, chunkSize);
}

template <class SignalR, class SignalU, class SignalV, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalU, SignalV>, int> = 0>
auto Filter(SignalR&& out, const SignalU& signal, const SignalV& filter, impl::ConvCentral, impl::FilterConv) {
	Convolution(out, signal, filter, CONV_CENTRAL);
}

template <class SignalR, class SignalU, class SignalV, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalU, SignalV>, int> = 0>
auto Filter(SignalR&& out, const SignalU& signal, const SignalV& filter, impl::ConvFull, impl::FilterOla, size_t chunkSize = 0) {
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
auto Filter(SignalR&& out, const SignalU& signal, const SignalV& filter, SignalS& state, impl::FilterOla, size_t chunkSize = 0) {
	assert(state.size() == filter.size() - 1);
	assert(out.size() == signal.size());

	std::fill(out.begin(), out.end(), remove_complex_t<typename std::decay_t<SignalR>::value_type>(0));
	OverlapAdd(AsView(out).subsignal(0, std::min(out.size(), state.size())), state, filter, filter.size() - 1, chunkSize, false);
	OverlapAdd(out, signal, filter, 0, chunkSize, false);
	impl::ShiftFilterState(state, signal);
}

template <class SignalR,
		  class SignalU,
		  class SignalV,
		  class SignalS,
		  std::enable_if_t<is_mutable_signal_v<SignalR> && is_mutable_signal_v<SignalS> && is_same_domain_v<SignalR, SignalU, SignalV, SignalS>, int> = 0>
auto Filter(SignalR&& out, const SignalU& signal, const SignalV& filter, SignalS& state, impl::FilterConv) {
	assert(state.size() == filter.size() - 1);
	assert(out.size() == signal.size());

	std::fill(out.begin(), out.end(), remove_complex_t<typename std::decay_t<SignalR>::value_type>(0));
	Convolution(AsView(out).subsignal(0, std::min(out.size(), state.size())), state, filter, filter.size() - 1, false);
	Convolution(out, signal, filter, 0, false);
	impl::ShiftFilterState(state, signal);
}

template <class SignalU, class SignalV, std::enable_if_t<is_same_domain_v<SignalU, SignalV>, int> = 0>
auto Filter(const SignalU& signal, const SignalV& filter, impl::ConvCentral, impl::FilterOla, size_t chunkSize = 0) {
	impl::ProductSignalT<SignalU, SignalV> out(ConvolutionLength(signal.size(), filter.size(), CONV_CENTRAL));
	Filter(out, signal, filter, CONV_CENTRAL, FILTER_OLA, chunkSize);
	return out;
}

template <class SignalU, class SignalV, std::enable_if_t<is_same_domain_v<SignalU, SignalV>, int> = 0>
auto Filter(const SignalU& signal, const SignalV& filter, impl::ConvCentral, impl::FilterConv) {
	impl::ProductSignalT<SignalU, SignalV> out(ConvolutionLength(signal.size(), filter.size(), CONV_CENTRAL));
	Filter(out, signal, filter, CONV_CENTRAL, FILTER_CONV);
	return out;
}

template <class SignalU, class SignalV, std::enable_if_t<is_same_domain_v<SignalU, SignalV>, int> = 0>
auto Filter(const SignalU& signal, const SignalV& filter, impl::ConvFull, impl::FilterOla, size_t chunkSize = 0) {
	impl::ProductSignalT<SignalU, SignalV> out(ConvolutionLength(signal.size(), filter.size(), CONV_FULL));
	Filter(out, signal, filter, CONV_FULL, FILTER_OLA, chunkSize);
	return out;
}

template <class SignalU, class SignalV, std::enable_if_t<is_same_domain_v<SignalU, SignalV>, int> = 0>
auto Filter(const SignalU& signal, const SignalV& filter, impl::ConvFull, impl::FilterConv) {
	impl::ProductSignalT<SignalU, SignalV> out(ConvolutionLength(signal.size(), filter.size(), CONV_FULL));
	Filter(out, signal, filter, CONV_FULL, FILTER_CONV);
	return out;
}

template <class SignalU,
		  class SignalV,
		  class SignalS,
		  std::enable_if_t<is_mutable_signal_v<SignalS> && is_same_domain_v<SignalU, SignalV, SignalS>, int> = 0>
auto Filter(const SignalU& signal, const SignalV& filter, SignalS&& state, impl::FilterOla, size_t chunkSize = 0) {
	impl::ProductSignalT<SignalU, SignalV> out(signal.size());
	Filter(out, signal, filter, state, FILTER_OLA, chunkSize);
	return out;
}

template <class SignalU,
		  class SignalV,
		  class SignalS,
		  std::enable_if_t<is_mutable_signal_v<SignalS> && is_same_domain_v<SignalU, SignalV, SignalS>, int> = 0>
auto Filter(const SignalU& signal, const SignalV& filter, SignalS&& state, impl::FilterConv) {
	impl::ProductSignalT<SignalU, SignalV> out(signal.size());
	Filter(out, signal, filter, state, FILTER_CONV);
	return out;
}

} // namespace dspbb