#pragma once

#include "../../LTISystems/TransferFunctionSystem.hpp"
#include "Realizations.hpp"


namespace dspbb {

template <class SignalT, class T, class U>
auto Filter(const SignalT& signal, const DiscreteTransferFunctionSystem<U>& filter, DirectFormI<T>& state) {
	SignalT r;
	r.Reserve(signal.Size());
	for (auto& sample : signal) {
		r.PushBack(state.Feed(sample, filter.Numerator(), filter.Denominator()));
	}
	return r;
}

} // namespace dspbb