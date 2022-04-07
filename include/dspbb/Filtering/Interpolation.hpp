#pragma once

#include "../Math/DotProduct.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalTraits.hpp"
#include "../Primitives/SignalView.hpp"
#include "Polyphase.hpp"

#include <numeric>


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
	using T = std::remove_const_t<typename signal_traits<SignalT>::type>;
	constexpr auto domain = signal_traits<SignalT>::domain;
	BasicSignal<T, domain> output((input.Size() + rate - 1) / rate);
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
	using T = std::remove_const_t<typename signal_traits<SignalT>::type>;
	constexpr auto domain = signal_traits<SignalT>::domain;
	BasicSignal<T, domain> output(input.Size() * rate);
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
		  class P,
		  eSignalDomain D,
		  std::enable_if_t<is_same_domain_v<SignalR, SignalT, BasicSignal<P, D>> && is_mutable_signal_v<SignalR>, int> = 0>
void Interpolate(SignalR&& hrOutput,
				 const SignalT& lrInput,
				 const PolyphaseView<P, D>& polyphase,
				 size_t hrOffset) {
	const ptrdiff_t rate = polyphase.FilterCount();
	const ptrdiff_t hrFilterSize = polyphase.OriginalSize();
	const ptrdiff_t lrPhaseSize = polyphase.PhaseSize();
	const ptrdiff_t hrOutputSize = hrOutput.Size();
	const ptrdiff_t lrInputSize = lrInput.Size();
	const ptrdiff_t hrInputSize = lrInputSize * rate;

	const ptrdiff_t hrOutputMaxSize = ConvolutionLength(hrInputSize, hrFilterSize, CONV_FULL);
	assert(hrOffset + hrOutputSize <= hrOutputMaxSize);

	for (size_t hrOutputIdx = hrOffset; hrOutputIdx < hrOffset + hrOutputSize; ++hrOutputIdx) {
		const ptrdiff_t hrInputIdx = 1 - hrFilterSize + hrOutputIdx;
		const ptrdiff_t lrInputIdx = (hrInputIdx + hrFilterSize - 1) / rate - lrPhaseSize + 1;
		const ptrdiff_t polyphaseIdx = (hrInputIdx + hrFilterSize - 1) % rate;

		const auto& phase = polyphase[polyphaseIdx];

		const Interval inputSpan = { ptrdiff_t(0), ptrdiff_t(lrInput.Size()) };
		const Interval lrInputInterval = { lrInputIdx, lrInputIdx + lrPhaseSize };
		const Interval lrPhaseInterval = { lrInputInterval.last - ptrdiff_t(phase.Size()), lrInputInterval.last };
		const Interval lrInputProductInterval = Intersection(inputSpan, Intersection(lrInputInterval, lrPhaseInterval));
		const Interval lrPhaseProductInterval = lrInputProductInterval - lrInputIdx;

		if (lrInputProductInterval.Size() > 0) {
			const auto lrInputView = AsView(lrInput).SubSignal(lrInputProductInterval.first,
															   lrInputProductInterval.last - lrInputProductInterval.first);
			const auto lrPhaseView = phase.SubSignal(lrPhaseProductInterval.first - lrPhaseSize + ptrdiff_t(phase.Size()),
													 lrPhaseProductInterval.last - lrPhaseProductInterval.first);
			const auto value = DotProduct(lrInputView, lrPhaseView);
			hrOutput[hrOutputIdx - hrOffset] = value;
		}
	}
}

template <class SignalT, class P, eSignalDomain Domain, std::enable_if_t<is_same_domain_v<SignalT, BasicSignal<P, Domain>>, int> = 0>
auto Interpolate(const SignalT& lrInput,
				 const PolyphaseView<P, Domain>& polyphase,
				 size_t hrOffset,
				 size_t hrLength) {
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	using R = multiplies_result_t<T, P>;

	BasicSignal<R, Domain> out(hrLength, R(0));
	Interpolate(out, lrInput, polyphase, hrOffset);
	return out;
}


namespace impl {
	template <class T1, class T2>
	constexpr auto lcm(T1 x, T2 y) {
		return std::lcm(x, y);
	}
	template <class T1, class... T>
	constexpr auto lcm(T1 head, T... tail) {
		return std::lcm(head, dspbb::impl::lcm(tail...));
	}
} // namespace impl


namespace resample {

	template <class ConvType>
	constexpr std::pair<uint64_t, uint64_t> ResamplingLength(size_t inputSize,
															 size_t filterSize,
															 size_t numPhases,
															 std::pair<uint64_t, uint64_t> sampleRates,
															 const ConvType&) {
		static_assert(std::is_same_v<ConvType, impl::ConvFull> || std::is_same_v<ConvType, impl::ConvCentral>);
		const uint64_t interpolatedSize = numPhases * inputSize;
		const uint64_t filteredSize = ConvolutionLength(interpolatedSize, filterSize, ConvType{});

		const uint64_t sampleRatesGcd = std::gcd(sampleRates.first, sampleRates.second);
		const uint64_t sampleCountGcd = std::gcd(filteredSize, sampleRates.first / sampleRatesGcd);
		const std::pair length = {
			(filteredSize / sampleCountGcd) * (sampleRates.second / sampleRatesGcd),
			sampleRates.first * numPhases / sampleCountGcd / sampleRatesGcd
		};
		return length;
	}

	template <class ConvType>
	constexpr std::pair<uint64_t, uint64_t> ResamplingStartPoint(size_t filterSize,
																 size_t numPhases,
																 std::pair<uint64_t, uint64_t> sampleRates,
																 const ConvType&) {
		if constexpr (std::is_same_v<ConvType, impl::ConvFull>) {
			return { 0, 1 };
		}
		else if constexpr (std::is_same_v<ConvType, impl::ConvFull>) {
			throw std::logic_error("not implemented");
		}
		else {
			static_assert(false && sizeof(ConvType), "Please use either CONV_FULL or CONV_CENTRAL");
		}
	}

	constexpr std::pair<uint64_t, uint64_t> ChangeSampleRate(std::pair<uint64_t, uint64_t> sampleRates,
															 std::pair<uint64_t, uint64_t> sample,
															 bool simplify = false) {
		if (simplify) {
			const auto gcdOutput = std::gcd(sample.first, sample.second);
			const auto outputIndexSimplified = std::pair{ sample.first / gcdOutput, sample.second / gcdOutput };
			const auto gcdRate = std::gcd(sampleRates.second, sampleRates.first);
			const auto sampleRateSimplified = std::pair{ sampleRates.first / gcdRate, sampleRates.second / gcdRate };
			const auto gcd1 = std::gcd(sampleRateSimplified.second, outputIndexSimplified.second);
			const auto gcd2 = std::gcd(sampleRateSimplified.first, outputIndexSimplified.first);
			return {
				(outputIndexSimplified.first / gcd2) * (sampleRateSimplified.second / gcd1),
				(outputIndexSimplified.second / gcd1) * (sampleRateSimplified.first / gcd2)
			};
		}
		return {
			sample.first * sampleRates.second,
			sample.second * sampleRates.first
		};
	}

	struct PhaseSample {
		size_t inputIndex;
		size_t phaseIndex;
		uint64_t weight;
	};

	constexpr std::pair<PhaseSample, PhaseSample> InputIndex2Sample(std::pair<uint64_t, uint64_t> inputIndex, size_t numPhases) {
		const uint64_t firstIndex = inputIndex.first / inputIndex.second;
		const std::pair fractionalIndex = {
			inputIndex.first % inputIndex.second,
			inputIndex.second
		};
		const size_t firstPhase = fractionalIndex.first * numPhases / fractionalIndex.second;
		const size_t secondWeight = fractionalIndex.first * numPhases % fractionalIndex.second;
		const size_t secondPhase = (firstPhase + 1) % numPhases;
		const size_t firstWeight = fractionalIndex.second - secondWeight;
		const size_t secondIndex = secondPhase == 0 ? firstIndex + 1 : firstIndex;

		return {
			PhaseSample{ firstIndex, firstPhase, firstWeight },
			PhaseSample{ secondIndex, secondPhase, secondWeight }
		};
	}

	template <class SignalT, class SignalU>
	auto DotProductSample(const SignalT& input, const SignalU& filter, size_t inputReverseFirst) {
		const ptrdiff_t desiredFirst = ptrdiff_t(inputReverseFirst) - filter.Size() + 1;
		const ptrdiff_t desiredLast = ptrdiff_t(inputReverseFirst) + 1;
		const ptrdiff_t possibleFirst = std::max(ptrdiff_t(0), desiredFirst);
		const ptrdiff_t possibleLast = std::min(ptrdiff_t(input.Size()), desiredLast);
		const ptrdiff_t count = possibleLast - possibleFirst;
		assert(count >= 0);
		const ptrdiff_t offset = possibleFirst - desiredFirst;

		const auto inputView = AsConstView(input).SubSignal(possibleFirst, count);
		const auto filterView = AsConstView(filter).SubSignal(offset, count);
		return DotProduct(inputView, filterView);
	}

	constexpr double ResamplingFilterCutoff(std::pair<uint64_t, uint64_t> sampleRates, size_t numPhases) {
		const double base = 1.0 / double(numPhases);
		const double rate = std::min(1.0, double(sampleRates.second) / double(sampleRates.first));
		return base * rate;
	}

	constexpr std::pair<uint64_t, uint64_t> ResamplingDelay(size_t filterSize,
															size_t numPhases,
															std::pair<uint64_t, uint64_t> sampleRates) {
		return {
			(filterSize - 1) * sampleRates.second,
			sampleRates.first * 2 * numPhases
		};
	}

	struct ContinuationParams {
		size_t firstInputSample;
		std::pair<uint64_t, uint64_t> startPoint;
	};

	constexpr ContinuationParams Continuation(std::pair<uint64_t, uint64_t> nextOutputSample,
											  size_t filterSize,
											  size_t numPhases,
											  std::pair<uint64_t, uint64_t> sampleRates) {
		const auto nextInputSample = ChangeSampleRate({ sampleRates.second, sampleRates.first }, nextOutputSample, true);
		const auto nextInterpolatedSample = std::pair{ nextInputSample.first * numPhases, nextInputSample.second * numPhases };
		const auto firstInputSample = std::pair{ int64_t(nextInterpolatedSample.first) - int64_t(nextInputSample.second * filterSize) + int64_t(nextInputSample.second),
												 nextInterpolatedSample.second };
		if (firstInputSample.first <= 0) {
			return { 0, nextOutputSample };
		}
		else {
			const size_t firstInputSampleWhole = firstInputSample.first / firstInputSample.second;
			const auto inputStartPoint = std::pair<uint64_t, uint64_t>{
				firstInputSample.first - firstInputSampleWhole * firstInputSample.second + int64_t(nextInputSample.second * filterSize) - int64_t(nextInputSample.second),
				firstInputSample.second
			};
			const auto outputStartPoint = ChangeSampleRate(sampleRates, inputStartPoint, true);
			return { firstInputSampleWhole, outputStartPoint };
		}
	}

} // namespace resample

using resample::ResamplingLength;
using resample::ResamplingFilterCutoff;
using resample::ResamplingDelay;


template <class SignalR,
		  class SignalT,
		  class P,
		  eSignalDomain D,
		  std::enable_if_t<is_same_domain_v<SignalR, SignalT, BasicSignal<P, D>> && is_mutable_signal_v<SignalR>, int> = 0>
resample::ContinuationParams Resample(SignalR&& output,
			  const SignalT& input,
			  const PolyphaseView<P, D>& polyphase,
			  std::pair<uint64_t, uint64_t> sampleRates,
			  std::pair<uint64_t, uint64_t> startPoint = { 0, 1 }) {
	assert(sampleRates.first > 0);
	assert(sampleRates.second > 0);
	assert(startPoint.second > 0);
	assert(polyphase.FilterCount() > 0);

	const auto maxLength = resample::ResamplingLength(input.Size(), polyphase.OriginalSize(), polyphase.FilterCount(), sampleRates, CONV_FULL);

	auto outputIndex = startPoint;
	for (auto outputIt = output.begin(); outputIt != output.end(); ++outputIt, outputIndex.first += outputIndex.second) {
		const auto inputIndex = resample::ChangeSampleRate({ sampleRates.second, sampleRates.first }, outputIndex);
		const auto [firstSampleLoc, secondSampleLoc] = resample::InputIndex2Sample(inputIndex, polyphase.FilterCount());
		const auto firstSampleVal = resample::DotProductSample(input, polyphase[firstSampleLoc.phaseIndex], firstSampleLoc.inputIndex);
		const auto secondSampleVal = resample::DotProductSample(input, polyphase[secondSampleLoc.phaseIndex], secondSampleLoc.inputIndex);
		using CommonType = decltype(firstSampleVal);
		*outputIt = (firstSampleVal * CommonType(firstSampleLoc.weight) + secondSampleVal * CommonType(secondSampleLoc.weight))
					/ (CommonType(firstSampleLoc.weight) + CommonType(secondSampleLoc.weight));
	}

	return resample::Continuation()
}

template <class SignalT,
		  class P,
		  eSignalDomain Domain,
		  std::enable_if_t<is_same_domain_v<SignalT, BasicSignal<P, Domain>>, int> = 0>
auto Resample(const SignalT& input,
			  const PolyphaseView<P, Domain>& polyphase,
			  std::pair<uint64_t, uint64_t> sampleRates,
			  std::pair<uint64_t, uint64_t> startPoint,
			  size_t orLength) {
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	using R = multiplies_result_t<T, P>;

	BasicSignal<R, Domain> out(orLength, R(0));
	Resample(out, input, polyphase, sampleRates, startPoint);
	return out;
}


} // namespace dspbb