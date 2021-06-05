#pragma once

#include "../Math/DotProduct.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "../Primitives/SignalTraits.hpp"
#include "PolyphaseFilter.hpp"


namespace dspbb {



template <class SignalR,
		  class SignalT,
		  std::enable_if_t<is_same_domain_v<SignalR, SignalT> && is_mutable_signal_v<SignalR>, int> = 0>
void InterpolateZeroFill(SignalR&& output,
						 const SignalT& input,
						 uint64_t rate) {
	assert(output.Size() == input.Size() * rate);
	auto writeIt = output.begin();
	for (auto& v : input) {
		*writeIt = v;
		++writeIt;
		for (uint64_t i = rate; i > 1; --i) {
			*writeIt = 0;
			++writeIt;
		}
	}
}


template <class SignalR,
		  class SignalT,
		  class P, eSignalDomain D,
		  std::enable_if_t<is_same_domain_v<SignalR, SignalT, Signal<P, D>> && is_mutable_signal_v<SignalR>, int> = 0>
std::pair<int64_t, uint64_t> Interpolate(SignalR&& output,
										 const SignalT& input,
										 const PolyphaseDecomposition<P, D>& polyphase,
										 std::pair<uint64_t, uint64_t> sampleRates,
										 std::pair<int64_t, uint64_t> startPoint = { 0, 1 }) {
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	const size_t commonRate = sampleRates.first * sampleRates.second * startPoint.second * polyphase.numFilters; // Use smallest common multiple instead.

	// All these are in common rate.
	const int64_t commonStartPoint = startPoint.first * commonRate / startPoint.second;
	const uint64_t commonStep = sampleRates.first * commonRate / sampleRates.second;
	int64_t commonSample = commonStartPoint;

	for (auto& o : output) {
		const int64_t commonFraction = commonSample % commonRate;

		const int64_t polyphaseIndex = commonFraction * polyphase.numFilters / commonRate;
		const int64_t polyphaseIndexNext = (polyphaseIndex + 1) % polyphase.numFilters;
		const int64_t polyphaseFraction = commonFraction * polyphase.numFilters % commonRate;

		const int64_t inputIndex = commonSample / commonRate - polyphase[polyphaseIndex].Size() + 1;
		const int64_t inputIndexNext = commonSample / commonRate - polyphase[polyphaseIndexNext].Size() + 1 + (polyphaseIndexNext == 0);

		int64_t inputIndexClamped = std::max(inputIndex, int64_t(0));
		int64_t inputIndexClampedNext = std::max(inputIndexNext, int64_t(0));

		const auto filter = polyphase[polyphaseIndex].SubSignal(inputIndexClamped - inputIndex);
		const auto filterNext = polyphase[polyphaseIndexNext].SubSignal(inputIndexClampedNext - inputIndexNext);
		const auto inputSection = AsView(input).SubSignal(inputIndexClamped);
		const auto inputSectionNext = AsView(input).SubSignal(inputIndexClampedNext);

		const auto sample = DotProduct(inputSection, filter, std::min(inputSection.Size(), filter.Size()));
		const auto sampleNext = DotProduct(inputSectionNext, filterNext, std::min(inputSectionNext.Size(), filterNext.Size()));
		o = (sample * R(commonRate - polyphaseFraction) + sampleNext * R(polyphaseFraction)) / R(commonRate);

		// Next sample
		commonSample += commonStep;
	}
	return { commonSample, commonRate };
}


} // namespace dspbb