#pragma once

#include "../../LTISystems/Systems.hpp"
#include "Realizations.hpp"


namespace dspbb {

namespace impl {

	template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT, class System, class State>
	auto Filter(SignalR&& out, const SignalT& signal, const System& filter, State& state) {
		assert(out.size() == signal.size());
		state.feed(signal.begin(), signal.end(), out.begin(), filter);
	}

} // namespace impl


/// <summary> Apply an IIR filter to a signal. </summary>
/// <param name="out"> Output buffer for the filtered signal. </param>
/// <param name="signal"> The signal to be filtered. </param>
/// <param name="filter"> The filter applied to the signal. </param>
/// <param name="state"> The state of the filter. </param>
template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT, class T, class U>
auto Filter(SignalR&& out, const SignalT& signal, const DiscreteTransferFunction<U>& filter, DirectFormI<T>& state) {
	impl::Filter(out, signal, filter, state);
}


/// <summary> Apply an IIR filter to a signal. </summary>
/// <param name="out"> Output buffer for the filtered signal. </param>
/// <param name="signal"> The signal to be filtered. </param>
/// <param name="filter"> The filter applied to the signal. </param>
/// <param name="state"> The state of the filter. </param>
template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT, class T, class U>
auto Filter(SignalR&& out, const SignalT& signal, const DiscreteTransferFunction<U>& filter, DirectFormII<T>& state) {
	impl::Filter(out, signal, filter, state);
}


/// <summary> Apply an IIR filter to a signal. </summary>
/// <param name="out"> Output buffer for the filtered signal. </param>
/// <param name="signal"> The signal to be filtered. </param>
/// <param name="filter"> The filter applied to the signal. </param>
/// <param name="state"> The state of the filter. </param>
template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT, class T, class U>
auto Filter(SignalR&& out, const SignalT& signal, const CascadedBiquad<U>& filter, CascadedForm<T>& state) {
	impl::Filter(out, signal, filter, state);
}


/// <summary> Apply an IIR filter to a signal. </summary>
/// <param name="signal"> The signal to be filtered. </param>
/// <param name="filter"> The filter applied to the signal. </param>
/// <param name="state"> The state of the filter. </param>
/// <returns> The filtered signal. </returns>
template <signal_or_view SignalT, class T, class U>
auto Filter(const SignalT& signal, const DiscreteTransferFunction<U>& filter, DirectFormI<T>& state) {
	SignalT out(signal.size());
	Filter(out, signal, filter, state);
	return out;
}


/// <summary> Apply an IIR filter to a signal. </summary>
/// <param name="signal"> The signal to be filtered. </param>
/// <param name="filter"> The filter applied to the signal. </param>
/// <param name="state"> The state of the filter. </param>
/// <returns> The filtered signal. </returns>
template <signal_or_view SignalT, class T, class U>
auto Filter(const SignalT& signal, const DiscreteTransferFunction<U>& filter, DirectFormII<T>& state) {
	SignalT out(signal.size());
	Filter(out, signal, filter, state);
	return out;
}


/// <summary> Apply an IIR filter to a signal. </summary>
/// <param name="signal"> The signal to be filtered. </param>
/// <param name="filter"> The filter applied to the signal. </param>
/// <param name="state"> The state of the filter. </param>
/// <returns> The filtered signal. </returns>
template <signal_or_view SignalT, class T, class U>
auto Filter(const SignalT& signal, const CascadedBiquad<U>& filter, CascadedForm<T>& state) {
	SignalT out(signal.size());
	Filter(out, signal, filter, state);
	return out;
}

} // namespace dspbb