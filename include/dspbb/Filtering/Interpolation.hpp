#pragma once

#include "../Math/DotProduct.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "PolyphaseFilter.hpp"


namespace dspbb {


/// <summary> Interpolates a signal using cubic polynomials. </summary>
/// <param name="input"> The signal you want to resample. </param>
/// <param name="sampleRateIn"> <paramref name="input"/>'s sample rate. </param>
/// <param name="sampleRateOut"> The result's sample rate. </param>
/// <param name="offset"> Where from the input stream to take the first output sample.
///		Zero means the first output sample is taken at exactly input[0], while a value of <paramref name="sampleRateOut">
///		means the sample is taken at exactly input[0]. Any positive value is valid. </param>
/// <returns> The interpolated samples, as many as fits in <paramref name="input"/>. </returns>
template <class T>
TimeSignal<T> InterpolateCubic(const TimeSignal<T>& input, const uint64_t sampleRateIn, const uint64_t sampleRateOut, const uint64_t offset) {
	TimeSignal<T> output;

	const uint64_t advance = sampleRateIn / sampleRateOut;
	const uint64_t carry = sampleRateIn % sampleRateOut;

	uint64_t inputIndexAccumulator = offset / sampleRateOut;
	uint64_t carryAccumulator = offset % sampleRateOut;

	while (inputIndexAccumulator + (carryAccumulator + sampleRateOut - 1) / sampleRateOut < input.Size()) {
		carryAccumulator += carry;
		inputIndexAccumulator += advance + carryAccumulator / sampleRateOut;
		carryAccumulator %= sampleRateOut;

		const T sample{}; //= impl::InterpolateSample(input, inputIndexAccumulator, float(carryAccumulator) / float(sampleRateOut));
		output.Append({ sample });
	}
	return output;
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

	const size_t commonRate = sampleRateIn * sampleRateOut * offsetDen; // Use GCD instead.
	auto Output2Common = [&](size_t outputIdx) { return outputIdx * commonRate / sampleRateOut; };
	auto Common2Input = [&](size_t commonIdx) { return std::pair{ commonIdx * sampleRateIn / commonRate, commonIdx * sampleRateIn % commonRate }; };

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