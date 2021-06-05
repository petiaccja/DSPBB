#pragma once

#include "../Math/DotProduct.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "PolyphaseFilter.hpp"


namespace dspbb {



template <class T, eSignalDomain Domain>
void InterpolateZeroFill(SignalView<T, Domain> output,
						 SignalView<const T, Domain> input,
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


template <class T, eSignalDomain Domain>
std::pair<int64_t, uint64_t> Interpolate(SignalView<T, Domain> output,
										 SignalView<const T, Domain> input,
										 const PolyphaseDecomposition<T, Domain>& polyphase,
										 std::pair<uint64_t, uint64_t> sampleRates,
										 std::pair<int64_t, uint64_t> startPoint = { 0, 1 }) {
	const size_t commonRate = sampleRates.first * sampleRates.second * startPoint.second * polyphase.numFilters; // Use smallest common multiple instead.

	// All these are in common rate.
	const int64_t commonStartPoint = startPoint.first * commonRate / startPoint.second;
	const uint64_t commonStep = sampleRates.first * commonRate / sampleRates.second;
	int64_t commonSample = commonStartPoint;
	size_t outputIdx = 0;

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
		const auto inputSection = input.SubSignal(inputIndexClamped);
		const auto inputSectionNext = input.SubSignal(inputIndexClampedNext);

		const Signal<T, Domain> filterDbg = { filter.begin(), filter.end() };
		const Signal<T, Domain> inputSectionDbg = { inputSection.begin(), inputSection.end() };
		const Signal<T, Domain> filterDbgNext = { filterNext.begin(), filterNext.end() };
		const Signal<T, Domain> inputSectionDbgNext = { inputSectionNext.begin(), inputSectionNext.end() };

		const T sample = DotProduct(inputSection, filter, std::min(inputSection.Size(), filter.Size()));
		const T sampleNext = DotProduct(inputSectionNext, filterNext, std::min(inputSectionNext.Size(), filterNext.Size()));
		o = (sample * T(commonRate - polyphaseFraction) + sampleNext * T(polyphaseFraction)) / T(commonRate);

		// Next sample
		commonSample += commonStep;
		++outputIdx;
	}
	return { commonSample, commonRate };
}


} // namespace dspbb