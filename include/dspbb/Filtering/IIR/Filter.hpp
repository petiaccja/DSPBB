#pragma once

#include "../../LTISystems/TransferFunctionSystem.hpp"
#include "Structures.hpp"


namespace dspbb {

template <class SignalT, class T, class U>
auto Filter(const SignalT& signal, const DiscreteTransferFunctionSystem<U>& filter, DirectForm1<T>& state) {
	SignalT r;
	r.Reserve(signal.Size());
	const auto normalization = *--filter.Denominator().end();
	for (auto& sample : signal) {
		state.zeroState[0] = sample;
		std::rotate(state.zeroState.begin(), ++state.zeroState.begin(), state.zeroState.end());
		
		const auto numerator = SignalView<const U, eSignalDomain::DOMAINLESS>{ filter.Numerator().begin(), filter.Numerator().end() };
		const auto denominator = SignalView<const U, eSignalDomain::DOMAINLESS>{ filter.Denominator().begin(), --filter.Denominator().end() };
		const float forward = DotProduct(numerator, state.zeroState);
		const float recursive = DotProduct(denominator, state.poleState);
		const float out = (forward - recursive) / normalization;

		if (state.poleState.Size() > 0) {
			state.poleState[0] = out;
			std::rotate(state.poleState.begin(), ++state.poleState.begin(), state.poleState.end());
		}
		r.PushBack(out);
	}
	return r;
}

} // namespace dspbb