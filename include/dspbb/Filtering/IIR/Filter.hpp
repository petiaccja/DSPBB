#pragma once

#include "../../LTISystems/System.hpp"
#include "Realizations.hpp"


namespace dspbb {

template <class SignalT, class T, class U>
auto Filter(const SignalT& signal, const DiscreteTransferFunction<U>& filter, DirectFormI<T>& state) {
	SignalT r;
	r.Reserve(signal.Size());
	for (auto& sample : signal) {
		r.PushBack(state.Feed(sample, filter.numerator.Coefficients(), filter.denominator.Coefficients()));
	}
	return r;
}

} // namespace dspbb