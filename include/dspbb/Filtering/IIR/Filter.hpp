#pragma once

#include "../../LTISystems/Systems.hpp"
#include "Realizations.hpp"


namespace dspbb {

namespace impl {
	template <class SignalR, class SignalT, class System, class State, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
	auto Filter(SignalR&& out, const SignalT& signal, const System& filter, State& state) {
		if (out.Size() != signal.Size()) {
			throw std::invalid_argument("Output and input signals must have the same size.");
		}
		auto signalIt = signal.begin();
		auto outIt = out.begin();
		for (; outIt != out.end() && signalIt != signal.end(); ++outIt, ++signalIt) {
			*outIt = state.Feed(*signalIt, filter);
		}
	}
} // namespace impl

template <class SignalR, class SignalT, class T, class U, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
auto Filter(SignalR&& out, const SignalT& signal, const DiscreteTransferFunction<U>& filter, DirectFormI<T>& state) {
	impl::Filter(out, signal, filter, state);
}

template <class SignalR, class SignalT, class T, class U, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
auto Filter(SignalR&& out, const SignalT& signal, const DiscreteTransferFunction<U>& filter, DirectFormII<T>& state) {
	impl::Filter(out, signal, filter, state);
}

template <class SignalR, class SignalT, class T, class U, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
auto Filter(SignalR&& out, const SignalT& signal, const CascadedBiquad<U>& filter, CascadedForm<T>& state) {
	impl::Filter(out, signal, filter, state);
}

template <class SignalT, class T, class U>
auto Filter(const SignalT& signal, const DiscreteTransferFunction<U>& filter, DirectFormI<T>& state) {
	SignalT out(signal.Size());
	Filter(out, signal, filter, state);
	return out;
}

template <class SignalT, class T, class U>
auto Filter(const SignalT& signal, const DiscreteTransferFunction<U>& filter, DirectFormII<T>& state) {
	SignalT out(signal.Size());
	Filter(out, signal, filter, state);
	return out;
}

template <class SignalT, class T, class U>
auto Filter(const SignalT& signal, const CascadedBiquad<U>& filter, CascadedForm<T>& state) {
	SignalT out(signal.Size());
	Filter(out, signal, filter, state);
	return out;
}

} // namespace dspbb