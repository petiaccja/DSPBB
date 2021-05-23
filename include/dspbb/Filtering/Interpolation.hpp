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


/// <summary> Interpolates a signal using cubic polynomials. </summary>
/// <param name="input"> The signal you want to resample. </param>
/// <param name="filter"> The polyphase filter to use for interpolation. See remarks. </param>
/// <param name="sampleRateIn"> <paramref name="input"/>'s sample rate. </param>
/// <param name="sampleRateOut"> The result's sample rate. </param>
/// <param name="offset"> Where from the input stream to take the first output sample.
///		Zero means the first output sample is taken at exactly input[0], while a value of <paramref name="sampleRateOut">
///		means the sample is taken at exactly input[1]. Any positive value is valid. </param>
/// <returns> The interpolated samples, as many as fits in <paramref name="input"/>. </returns>
/// <remarks> The polyphase filter can have an arbitrary number of filters. It's cutoff frequency
///		must be at most the lower of the Nyquist frequencies of either the input or output
///		sample rates. </remarks>
template <class T>
TimeSignal<T> InterpolatePolyphase(SignalView<const T, TIME_DOMAIN> input,
								   const PolyphaseFilter<T>& filter,
								   const uint64_t sampleRateIn,
								   const uint64_t sampleRateOut,
								   const uint64_t offsetNum,
								   const uint64_t offsetDen) {

	TimeSignal<T> output;
	const size_t polyphaseCount = filter.NumFilters();
	const size_t polyphaseTaps = filter.NumTaps();

	//const size_t commonRate = sampleRateIn * sampleRateOut * offsetDen; // Use GCD instead.
	//auto Output2Common = [&](size_t outputIdx) { return outputIdx * commonRate / sampleRateOut; };
	//auto Common2Input = [&](size_t commonIdx) { return std::pair<uint64_t, uint64_t>{ commonIdx * sampleRateIn / commonRate, commonIdx * sampleRateIn % commonRate }; };

	size_t outputIdx = 0;
	while (true) {
		const uint64_t inputIdxHiRate = outputIdx * sampleRateIn;
		const uint64_t inputIdx = inputIdxHiRate / sampleRateOut;
		const uint64_t inputIdxFractionHiRate = inputIdxHiRate % sampleRateOut;

		const uint64_t leftPolyphaseIdx = inputIdxFractionHiRate * polyphaseCount / sampleRateOut;
		const uint64_t polyphaseFractionHiRate = inputIdxFractionHiRate * polyphaseCount % sampleRateOut;
		const uint64_t rightPolyphaseIdx = (leftPolyphaseIdx + 1) % polyphaseCount;
		const uint64_t rightInputIdx = inputIdx + size_t(rightPolyphaseIdx == 0);

		if (rightInputIdx + polyphaseTaps >= input.Size()) {
			break;
		}

		const T sampleLeft = DotProduct(input.SubSignal(inputIdx), filter.Filter(leftPolyphaseIdx), polyphaseTaps);
		const T sampleRight = DotProduct(input.SubSignal(rightInputIdx), filter.Filter(rightPolyphaseIdx), polyphaseTaps);
		const T lerpParam = T(polyphaseFractionHiRate) / T(sampleRateOut);
		const T sample = sampleLeft * (T(1) - lerpParam) + sampleRight * lerpParam;

		output.Append({ sample });
		++outputIdx;
	}
	return output;
}


} // namespace dspbb