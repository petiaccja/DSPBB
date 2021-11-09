#pragma once

#include "../../LTISystems/Systems.hpp"
#include "Realizations.hpp"


namespace dspbb {

template <class SignalT, class T, class U>
auto Filter(const SignalT& signal, const DiscreteTransferFunction<U>& filter, DirectFormI<T>& state) {
	SignalT r;
	r.Reserve(signal.Size());
	for (auto& sample : signal) {
		r.PushBack(state.Feed(sample, filter));
	}
	return r;
}

template <class SignalT, class T, class U>
auto Filter(const SignalT& signal, const DiscreteTransferFunction<U>& filter, DirectFormII<T>& state) {
	SignalT r;
	r.Reserve(signal.Size());
	for (auto& sample : signal) {
		r.PushBack(state.Feed(sample, filter));
	}
	return r;
}

template <class SignalT, class T, class U>
auto Filter(const SignalT& signal, const CascadedBiquad<U>& filter, CascadedForm<T>& state) {
	SignalT r;
	r.Reserve(signal.Size());
	for (auto& sample : signal) {
		r.PushBack(state.Feed(sample, filter));
	}
	return r;
}

} // namespace dspbb