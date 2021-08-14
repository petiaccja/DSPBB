#pragma once

#include "../Math/DotProduct.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalTraits.hpp"
#include "../Primitives/SignalView.hpp"
#include "PolyphaseFilter.hpp"


namespace dspbb {


/// <summary>
/// Erases all but every <paramref name="rate"/>th sample.
/// </summary>
template <class SignalR,
		  class SignalT,
		  std::enable_if_t<is_same_domain_v<SignalR, SignalT> && is_mutable_signal_v<SignalR>, int> = 0>
void Decimate(SignalR&& output,
			  const SignalT& input,
			  size_t rate) {
	assert(output.Size() == (input.Size() + rate - 1) / rate);
	size_t readIdx = 0;
	for (auto& o : output) {
		o = input[readIdx];
		readIdx += rate;
	}
}

template <class SignalT, std::enable_if_t<is_signal_like_v<SignalT>, int> = 0>
auto Decimate(const SignalT& input, size_t rate) {
	SignalT output((input.Size() + rate - 1) / rate);
	Decimate(output, input, rate);
	return output;
}


/// <summary>
/// Inserts zeros between samples to increase sample rate by a factor of <paramref name="rate"/>.
/// </summary>
/// <remarks> Follow expansion by a low-pass filter to interpolate a signal. </remarks>
template <class SignalR,
		  class SignalT,
		  std::enable_if_t<is_same_domain_v<SignalR, SignalT> && is_mutable_signal_v<SignalR>, int> = 0>
void Expand(SignalR&& output,
			const SignalT& input,
			size_t rate) {
	assert(output.Size() == input.Size() * rate);
	auto writeIt = output.begin();
	for (auto& i : input) {
		*writeIt = i;
		++writeIt;
		for (size_t i = rate; i > 1; --i) {
			*writeIt = 0;
			++writeIt;
		}
	}
}

template <class SignalT, std::enable_if_t<is_signal_like_v<SignalT>, int> = 0>
auto Expand(const SignalT& input, size_t rate) {
	SignalT output(input.Size() * rate);
	Expand(output, input, rate);
	return output;
}


/// <summary>
/// Inserts meaningful samples to increase sample rate by a factor of <paramref name="polyphase"/>.numFilters.
/// </summary>
/// <param name="polyphase"> A polyphase decomposition of an appropriate low-pass filter.
///		The number of phases defines the interpolation ratio. </param>
///	<remarks> No need to follow up with low-pass filtering.
///		The polyphase filter must have the appropriate cutoff-frequency of
///		(input sample rate / 2), and the polyphase filter must operate at
///		the output sample rate. </remarks>
template <class SignalR,
		  class SignalT,
		  class P, eSignalDomain D,
		  std::enable_if_t<is_same_domain_v<SignalR, SignalT, Signal<P, D>> && is_mutable_signal_v<SignalR>, int> = 0>
void Interpolate(SignalR&& output,
				 const SignalT& input,
				 const PolyphaseDecomposition<P, D>& polyphase,
				 intptr_t offset) {
	const auto rate = polyphase.numFilters;
	size_t count = output.Size() / rate;

	auto writeIt = output.begin();
	for (size_t inputIdx = 0; inputIdx < count; ++inputIdx) {
		for (size_t i = 0; i < rate; ++i, ++writeIt) {
			const auto& filter = polyphase[i];
			intptr_t offsetInputIdx = intptr_t(inputIdx) - intptr_t(filter.Size()) + 1 + intptr_t(offset);
			const intptr_t inputIndexClamped = std::max(offsetInputIdx, intptr_t(0));
			const auto subInput = AsConstView(input).SubSignal(inputIndexClamped);
			const auto subFilter = filter.SubSignal(-std::min(intptr_t(0), offsetInputIdx));
			const size_t dotLength = std::min(subFilter.Size(), subInput.Size());
			*writeIt = DotProduct(subFilter.SubSignal(0, dotLength), subInput.SubSignal(0, dotLength));
		}
	}
}


/// <summary>
/// Arbitrary rational resampling of a signal using approximate polyphase interpolation.
/// </summary>
/// <param name="polyphase"> A polyphase decomposition of an appropriate low-pass filter.
///		The number of phases is arbitrary (see remarks), but must be at least 2. </param>
/// <param name="sampleRates"> A pair of the {input, output} sample rates.
///		For example, use {16000, 44100} if you want to increase the sample rate of a signal by 2.75625 times. </param>
/// <param name="startPoint"> A rational number that tells where to take the first output sample. For example, specify {1, 2}
///		so that the first output sample will be taken from halfway between the first and second input samples. Specify {3, 77}
///		to take the first output sample from just after the first input sample. The denominator is arbitrary, the numerator
///		must be larger than zero, and float(num)/float(den) must be smaller or equal to length(<paramref name="input"/>)-1.</param>
/// <returns> A rational number that tells from where the next output sample would have been taken.
///		Analogues in meaning to <paramref name="startPoint"/>, but its denominator is not necessarily the same.
///		When resampling streams in batches, you can use this as a base for the <paramref name="startPoint"/> to a subsequent Resample call to ensure
///		the output stream is continuous. </returns>
template <class SignalR,
		  class SignalT,
		  class P, eSignalDomain D,
		  std::enable_if_t<is_same_domain_v<SignalR, SignalT, Signal<P, D>> && is_mutable_signal_v<SignalR>, int> = 0>
std::pair<int64_t, uint64_t> Resample(SignalR&& output,
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

		const auto dotLength = std::min(inputSection.Size(), filter.Size());
		const auto sample = DotProduct(inputSection.SubSignal(0, dotLength), filter.SubSignal(0, dotLength));
		const auto dotLengthNext = std::min(inputSectionNext.Size(), filterNext.Size());
		const auto sampleNext = DotProduct(inputSectionNext.SubSignal(0, dotLengthNext), filterNext.SubSignal(0, dotLengthNext));
		o = (sample * R(commonRate - polyphaseFraction) + sampleNext * R(polyphaseFraction)) / R(commonRate);

		// Next sample
		commonSample += commonStep;
	}
	return { commonSample, commonRate };
}


} // namespace dspbb