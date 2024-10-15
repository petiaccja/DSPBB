#pragma once

#include "../Generators/Spaces.hpp"
#include "../Math/DotProduct.hpp"
#include "../Math/Rational.hpp"
#include "../Signal/Signal.hpp"
#include "../Signal/SignalView.hpp"
#include "../Signal/Traits.hpp"
#include "Polyphase.hpp"

namespace dspbb {

//------------------------------------------------------------------------------
// Public utilities
//------------------------------------------------------------------------------

/// <summary> The parameters necessary for interpolating the next chunk of the signal. </summary>
struct InterpolSuspensionPoint {
	/// <summary> This is the first sample from the currently processed input signal
	///		that needs to be prepended to the next chunk of input signal. </summary>
	size_t firstInputSample;
	/// <summary> This is the value that you have to pass for offset to the interpolation
	///		of the next chunk. </summary>
	size_t outputOffset;
};


/// <summary> The parameters necessary for resampling the next chunk of the signal. </summary>
struct ResampleSuspensionPoint {
	/// <summary> This is the first sample from the currently processed input signal
	///		that needs to be prepended to the next chunk of input signal. </summary>
	size_t firstInputSample;
	/// <summary> This is the value that you have to pass for offset to the resampling
	///		of the next chunk. </summary>
	Rational<int64_t> outputOffset;
};


/// <summary> Compute the length of the interpolated signal. </summary>
/// <param name="inputSize"> The length of the signal to be resampled. </param>
/// <param name="filterSize"> The number of coefficients of the low-pass filter. </param>
/// <param name="factor"> The number of phases used for the polyphase decomposition -- the interpolation factor. </param>
/// <param name="convMethod"> Whether you want full or central convolution. </param>
constexpr size_t InterpolLength(size_t inputSize,
								size_t filterSize,
								size_t factor,
								eConvolutionMethod convMethod) {
	const ptrdiff_t expandedInputSize = inputSize * factor;
	return ConvolutionLength(expandedInputSize, filterSize, convMethod);
}


/// <summary> Compute the output offset for the interpolating functions for full and central convolution. </summary>
/// <param name="inputSize"> The length of the signal to be resampled. </param>
/// <param name="filterSize"> The number of coefficients of the low-pass filter. </param>
/// <param name="factor"> The number of phases used for the polyphase decomposition. </param>
/// <param name="convMethod"> Whether you want full or central convolution. </param>
constexpr size_t InterpolOffset(size_t inputSize,
								size_t filterSize,
								size_t factor,
								eConvolutionMethod convMethod) {
	return convMethod == eConvolutionMethod::CENTRAL ? (std::min(filterSize, inputSize * factor) - 1) : size_t(0);
}


/// <summary> Compute the normalized cutoff frequency of the low-pass filter for interpolation. </summary>
/// <param name="factor"> The interpolation factor. </param>
constexpr double InterpolFilterCutoff(size_t factor) {
	return 1.0 / double(factor);
}


/// <summary> Compute how many samples the output is delayed due to the filter. </summary>
/// <remarks>
/// <para> For simple convolution, this would normally be half the filter size
///		for a linear phase FIR filter. For interpolation, this function is more convenient. </para>
/// <para> The delay is expressed in output samples, at the output sample rate. </para>
///	</remarks>
constexpr size_t InterpolDelay(size_t filterSize) {
	return (filterSize - 1) / 2;
}


/// <summary> Compute the length of the resampled signal for full and central convolution. </summary>
/// <param name="inputSize"> The length of the signal to be resampled. </param>
/// <param name="filterSize"> The number of coefficients of the low-pass filter. </param>
/// <param name="numPhases"> The number of phases used for the polyphase decomposition. </param>
/// <param name="sampleRates"> The sample rates of the input and output signals. </param>
/// <param name="convMethod"> Whether you want full or central convolution. </param>
constexpr Rational<int64_t> ResampleLength(size_t inputSize,
										   size_t filterSize,
										   size_t numPhases,
										   Rational<int64_t> sampleRates,
										   eConvolutionMethod convMethod) {
	const auto expandedInputSize = int64_t(inputSize * numPhases);
	const auto filteredExpandedInputSize = int64_t(ConvolutionLength(expandedInputSize, filterSize, convMethod));

	return filteredExpandedInputSize / sampleRates / int64_t(numPhases);
}


/// <summary> Compute the output offset for the resampling functions for full and central convolution. </summary>
/// <param name="inputSize"> The length of the signal to be resampled. </param>
/// <param name="filterSize"> The number of coefficients of the low-pass filter. </param>
/// <param name="numPhases"> The number of phases used for the polyphase decomposition. </param>
/// <param name="sampleRates"> The sample rates of the input and output signals. </param>
/// <param name="convMethod"> Whether you want full or central convolution. </param>
constexpr Rational<int64_t> ResampleOffset(size_t inputSize,
										   size_t filterSize,
										   size_t numPhases,
										   Rational<int64_t> sampleRates,
										   eConvolutionMethod convMethod) {
	const Rational<int64_t> startPointIn = {
		int64_t(std::min(filterSize, inputSize * numPhases) - 1),
		int64_t(numPhases)
	};

	return convMethod == eConvolutionMethod::CENTRAL ? startPointIn / sampleRates : Rational<int64_t>{ 0 };
}


/// <summary> Compute the normalized cutoff frequency of the low-pass filter for resampling. </summary>
///	<param name="sampleRates"> The sample rate of the input and output signals. </param>
/// <param name="numPhases"> The number of the phases of the polyphase LPF used for the resampling. </param>
constexpr double ResampleFilterCutoff(Rational<int64_t> sampleRates, size_t numPhases) {
	const double base = 1.0 / double(numPhases);
	const double rate = std::min(1.0, 1.0 / double(sampleRates));
	return base * rate;
}


/// <summary> Compute how many samples the output is delayed due to the filter. </summary>
/// <remarks>
/// <para> For simple convolution, this would normally be half the filter size
///		for a linear phase FIR filter. For resampling, the formula is more complex. </para>
/// <para> The delay is expressed in output samples, at the output sample rate. </para>
///	</remarks>
constexpr Rational<int64_t> ResampleDelay(size_t filterSize,
										  size_t numPhases,
										  Rational<int64_t> sampleRates) {
	return Rational<int64_t>{ int64_t(filterSize) - 1, 2 * int64_t(numPhases) } / sampleRates;
}


//------------------------------------------------------------------------------
// Internal utilities
//------------------------------------------------------------------------------

namespace impl {
	inline InterpolSuspensionPoint FindInterpolSuspensionPoint(size_t nextOutputSample, size_t filterSize, size_t numPhases) {
		const ptrdiff_t firstOutputSample = ptrdiff_t(nextOutputSample) - ptrdiff_t(filterSize - 1);
		if (firstOutputSample < 0) {
			return { 0, nextOutputSample };
		}

		const size_t firstInputSample = firstOutputSample / numPhases;
		const size_t startPoint = firstOutputSample - (numPhases * firstInputSample) + (filterSize - 1);

		return { firstInputSample, startPoint };
	}

	constexpr Rational<int64_t> ChangeSampleRate(int64_t sourceRate,
												 int64_t targetRate,
												 Rational<int64_t> sample) {
		return sample * Rational{ targetRate, sourceRate };
	}

	constexpr ResampleSuspensionPoint FindResampleSuspensionPoint(Rational<int64_t> nextOutputSample,
																  size_t filterSize,
																  size_t numPhases,
																  Rational<int64_t> sampleRates) {
		const auto nextInputSample = ChangeSampleRate(sampleRates.Denominator(), sampleRates.Numerator(), nextOutputSample);
		const auto convolutionOffset = Rational{ int64_t(filterSize) - 1, int64_t(numPhases) };
		const auto firstInputSample = nextInputSample - convolutionOffset;

		if (firstInputSample <= 0) {
			return { 0, nextOutputSample };
		}
		else {
			const size_t firstInputSampleWhole = floor(firstInputSample);
			const auto inputStartPoint = frac(firstInputSample) + convolutionOffset;
			const auto outputStartPoint = ChangeSampleRate(sampleRates.Numerator(), sampleRates.Denominator(), inputStartPoint);
			return { firstInputSampleWhole, outputStartPoint };
		}
	}

	struct PhaseSample {
		size_t inputIndex;
		size_t phaseIndex;
		uint64_t weight;
	};

	constexpr std::pair<PhaseSample, PhaseSample> InputIndex2Sample(Rational<int64_t> inputIndex, size_t numPhases) {
		const Rational indexFrac = frac(inputIndex);

		const size_t firstPhase = floor(indexFrac * int64_t(numPhases));
		const size_t secondPhase = (firstPhase + 1) % numPhases;

		const Rational t = frac(indexFrac * int64_t(numPhases));
		const size_t secondWeight = t.Numerator();
		const size_t firstWeight = t.Denominator() - t.Numerator();

		const size_t firstIndex = floor(inputIndex);
		const size_t secondIndex = secondPhase == 0 ? firstIndex + 1 : firstIndex;

		return {
			PhaseSample{ firstIndex, firstPhase, firstWeight },
			PhaseSample{ secondIndex, secondPhase, secondWeight }
		};
	}

	template <class SignalT, class SignalU>
	auto DotProductSample(const SignalT& input, const SignalU& filter, size_t inputReverseFirst) {
		const ptrdiff_t desiredFirst = ptrdiff_t(inputReverseFirst) - filter.size() + 1;
		const ptrdiff_t desiredLast = ptrdiff_t(inputReverseFirst) + 1;
		const ptrdiff_t possibleFirst = std::max(ptrdiff_t(0), desiredFirst);
		const ptrdiff_t possibleLast = std::min(ptrdiff_t(input.size()), desiredLast);
		const ptrdiff_t count = possibleLast - possibleFirst;
		assert(count >= 0);
		const ptrdiff_t offset = possibleFirst - desiredFirst;

		const auto inputView = AsConstView(input).subsignal(possibleFirst, count);
		const auto filterView = AsConstView(filter).subsignal(offset, count);
		return DotProduct(inputView, filterView);
	}

} // namespace impl


//------------------------------------------------------------------------------
// Expansion & Interpolation & Resampling
//------------------------------------------------------------------------------


/// <summary> Keep only every <paramref name="factor"/>th sample. </summary>
/// <param name="output"> The decimated signal is written here. </param>
/// <param name="input"> The signal to decimate. </param>
/// <param name="factor"> The decimation factor. </param>
/// <remarks> This function does not do low-pass filtering. If you want to downsample,
///		you'll have to run a low-pass filter and then decimate. </remarks>
template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT>
void Decimate(SignalR&& output, const SignalT& input, size_t factor) {
	assert(!IsAliasing(output, input) || IsFullyAliasing(output, input));
	const auto count = std::min(output.size(), (input.size() + factor - 1) / factor);
	for (size_t i = 0; i < count; ++i) {
		output[i] = input[factor * i];
	}
}


/// <summary> Keep only every <paramref name="factor"/>th sample. </summary>
/// <param name="input"> The signal to decimate. </param>
/// <param name="factor"> The decimation factor. </param>
/// <remarks> This function does not low-pass filter the input. If you want to downsample,
///		you have to do low-pass filtering before decimation. </remarks>
template <signal_or_view SignalT>
auto Decimate(const SignalT& input, size_t factor) {
	using T = std::remove_const_t<scalar_type_t<SignalT>>;
	constexpr auto domain = domain_v<SignalT>;
	BasicSignal<T, domain> output((input.size() + factor - 1) / factor);
	Decimate(output, input, factor);
	return output;
}


/// <summary> Insert <paramref name="factor"/> - 1 zeros after every sample. </summary>
/// <param name="output"> The expanded signal is written here. </param>
/// <param name="input"> The signal to expand. </param>
/// <param name="factor"> The expansion factor. </param>
/// <remarks> This function does not low-pass filter the output. If you want to upsample,
///		you have to do low-pass filtering after expansion. </remarks>
template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT>
void Expand(SignalR&& output, const SignalT& input, size_t factor) {
	assert(!IsAliasing(output, input) || IsFullyAliasing(output, input));
	const auto count = std::min(output.size(), input.size() * factor);
	const auto zero = static_cast<scalar_type_t<std::decay_t<SignalR>>>(0);

	auto i = ptrdiff_t(count) - 1;
	auto j = ptrdiff_t((count + factor - 1) / factor) - 1;
	auto fraction = ptrdiff_t(i % factor);
	while (i >= 0) {
		output[i] = fraction == 0 ? input[j] : zero;
		i = i - 1;
		j = fraction == 0 ? j - 1 : j;
		fraction = fraction == 0 ? ptrdiff_t(factor) - 1 : fraction - 1;
	}
}


/// <summary> Insert <paramref name="factor"/> - 1 zeros after every sample. </summary>
/// <param name="input"> The signal to expand. </param>
/// <param name="factor"> The expansion factor. </param>
/// <remarks> This function does not low-pass filter the output. If you want to upsample,
///		you have to do low-pass filtering after expansion. </remarks>
template <signal_or_view SignalT>
auto Expand(const SignalT& input, size_t factor) {
	using T = std::remove_const_t<scalar_type_t<SignalT>>;
	constexpr auto domain = domain_v<SignalT>;
	BasicSignal<T, domain> output(input.size() * factor);
	Expand(output, input, factor);
	return output;
}


/// <summary> Interpolate (upsample) a signal using a polyphase FIR low-pass filter. </summary>
/// <param name="output"> The interpolated signal is written here. </param>
/// <param name="input"> The signal to interpolate. </param>
/// <param name="filter"> A polyphase low-pass filter. </param>
/// <param name="outputOffset"> Controls the starting point of the output. </param>
/// <returns> Information for interpolating the next chunk of input. </returns>
/// <remarks>
/// <para> The number of phases of <paramref name="filter"/> corresponds
///		to the interpolation factor. The filter should be a low-pass filter,
///		with a normalized cutoff frequency of 1/interpolation_factor. </para>
///	<para> Imagine the full (not central) interpolated output. The <paramref name="outputOffset"/>
///		gives the index to the first sample of this full output
///		that you want to write into <paramref name="output"/>. </para>
/// </remarks>
template <mutable_signal_or_view_r SignalR,
		  same_domain_as_r<SignalR> SignalT,
		  polyphase_or_view PolyphaseTy>
InterpolSuspensionPoint Interpolate(SignalR&& output,
									const SignalT& input,
									const PolyphaseTy& filter,
									size_t outputOffset) {
	assert(!IsAliasing(output, input));
	const ptrdiff_t rate = filter.num_phases();
	const ptrdiff_t hrFilterSize = filter.size_original();
	const ptrdiff_t lrPhaseSize = filter.size_per_phase();
	const ptrdiff_t hrOutputSize = output.size();

	const ptrdiff_t hrOutputMaxSize = InterpolLength(input.size(), hrFilterSize, rate, CONV_FULL);
	assert(ptrdiff_t(outputOffset) + hrOutputSize <= hrOutputMaxSize);

	size_t hrOutputIdx = outputOffset;
	for (; hrOutputIdx < outputOffset + hrOutputSize; ++hrOutputIdx) {
		const ptrdiff_t hrInputIdx = 1 - hrFilterSize + hrOutputIdx;
		const ptrdiff_t lrInputIdx = (hrInputIdx + hrFilterSize - 1) / rate - lrPhaseSize + 1;
		const ptrdiff_t polyphaseIdx = (hrInputIdx + hrFilterSize - 1) % rate;

		const auto& phase = filter[polyphaseIdx];

		const Interval inputSpan = { ptrdiff_t(0), ptrdiff_t(input.size()) };
		const Interval lrInputInterval = { lrInputIdx, lrInputIdx + lrPhaseSize };
		const Interval lrPhaseInterval = { lrInputInterval.last - ptrdiff_t(phase.size()), lrInputInterval.last };
		const Interval lrInputProductInterval = Intersection(inputSpan, Intersection(lrInputInterval, lrPhaseInterval));
		const Interval lrPhaseProductInterval = lrInputProductInterval - lrInputIdx;

		if (lrInputProductInterval.size() > 0) {
			const auto lrInputView = AsView(input).subsignal(lrInputProductInterval.first,
															 lrInputProductInterval.last - lrInputProductInterval.first);
			const auto lrPhaseView = phase.subsignal(lrPhaseProductInterval.first - lrPhaseSize + ptrdiff_t(phase.size()),
													 lrPhaseProductInterval.last - lrPhaseProductInterval.first);
			const auto value = DotProduct(lrInputView, lrPhaseView);
			output[hrOutputIdx - outputOffset] = value;
		}
	}

	return impl::FindInterpolSuspensionPoint(hrOutputIdx, filter.size_original(), filter.num_phases());
}


/// <summary> Interpolate (upsample) a signal using a polyphase FIR low-pass filter. </summary>
/// <param name="input"> The signal to interpolate. </param>
/// <param name="filter"> A polyphase low-pass filter. </param>
/// <param name="outputOffset"> Controls the starting point of the output. </param>
/// <param name="outputLength"> The number of interpolated output samples. </param>
/// <returns> The interpolated signal. </returns>
/// <remarks>
/// <para> The number of phases of <paramref name="filter"/> corresponds to the interpolation factor.
///		The low-pass cutoff frequency can be calculated with the utility function in this module. </para>
///	<para> Imagine the full (not central) interpolated output. The <paramref name="outputOffset"/>
///		gives the index to the first sample of this full output
///		that you want to write into <paramref name="output"/>. </para>
/// </remarks>
template <signal_or_view SignalTy, polyphase_or_view PolyphaseTy>
auto Interpolate(const SignalTy& input,
				 const PolyphaseTy& filter,
				 size_t outputOffset,
				 size_t outputLength) {
	using T = typename SignalTy::value_type;
	using U = typename PolyphaseTy::container::value_type;
	using R = multiplies_result_t<T, U>;
	constexpr auto Domain = domain_v<SignalTy>;

	BasicSignal<R, Domain> out(outputLength, R(0));
	const auto suspensionPoint = Interpolate(out, input, filter, outputOffset);
	return std::pair{ std::move(out), suspensionPoint };
}


/// <summary> Interpolate (upsample) a signal using a polyphase FIR low-pass filter. </summary>
/// <param name="output"> The interpolated signal is written here. </param>
/// <param name="input"> The signal to interpolate. </param>
/// <param name="filter"> The polyphase low-pass filter. </param>
/// <param name="convMethod"> Full or central convolution. </param>
/// <returns> The interpolated signal. </returns>
template <mutable_signal_or_view_r SignalOut, same_domain_as_r<SignalOut> SignalTy, polyphase_or_view PolyphaseTy, eConvolutionMethod ConvMethod>
auto Interpolate(const SignalOut&& output,
				 const SignalTy& input,
				 const PolyphaseTy& filter,
				 std::integral_constant<eConvolutionMethod, ConvMethod> convMethod) {
	const size_t outputLength = floor(InterpolLength(input.size(), filter.size_original(), filter.num_phases(), convMethod));
	const auto offset = InterpolOffset(input.size(), filter.size_original(), filter.num_phases(), convMethod);
	return Interpolate(output, input, filter, offset, outputLength);
}


/// <summary> Interpolate (upsample) a signal using a polyphase FIR low-pass filter. </summary>
/// <param name="input"> The signal to interpolate. </param>
/// <param name="filter"> The polyphase low-pass filter. </param>
/// <param name="convMethod"> Full or central convolution. </param>
/// <returns> The interpolated signal. </returns>
template <signal_or_view SignalTy, polyphase_or_view PolyphaseTy, eConvolutionMethod ConvMethod>
auto Interpolate(const SignalTy& input,
				 const PolyphaseTy& filter,
				 std::integral_constant<eConvolutionMethod, ConvMethod> convMethod) {
	const size_t outputLength = InterpolLength(input.size(), filter.size_original(), filter.num_phases(), convMethod);
	const auto offset = InterpolOffset(input.size(), filter.size_original(), filter.num_phases(), convMethod);
	return Interpolate(input, filter, offset, outputLength);
}


/// <summary> Resample the signal using a polyphase FIR low-pass filter. </summary>
/// <param name="output"> The resampled signal is written here. </param>
/// <param name="input"> The signal to resample. </param>
/// <param name="filter"> The polyphase low-pass filter. </param>
/// <param name="sampleRates"> The { input, output } sample rates. </param>
/// <param name="outputOffset"> Controls the starting point of the output. </param>
/// <returns> Information for resampling the next chunk of input. </returns>
/// <remarks>
/// <para> You're free to choose the polyphase filter's parameters, but the more phases
///		the smoother the output. The low-pass cutoff frequency can be calculated with the
///		utility function in this module. </para>
/// <para> Imagine the full (not central) resampled output. The <paramref name="outputOffset"/>
///		gives the fractional index to the first sample of this full output
///		that you want to write into <paramref name="output"/>. </para>
///	</remarks>
template <mutable_signal_or_view_r SignalR,
		  same_domain_as_r<SignalR> SignalT,
		  polyphase_or_view PolyphaseTy>
ResampleSuspensionPoint Resample(SignalR&& output,
								 const SignalT& input,
								 const PolyphaseTy& filter,
								 Rational<int64_t> sampleRates,
								 Rational<int64_t> outputOffset = { 0, 1 }) {
	assert(!IsAliasing(output, input));
	assert(sampleRates >= 0ll);
	assert(outputOffset >= 0ll);
	assert(filter.num_phases() > 0);

	[[maybe_unused]] const auto maxLength = ResampleLength(input.size(), filter.size_original(), filter.num_phases(), sampleRates, CONV_FULL);
	assert(outputOffset + int64_t(output.size()) <= maxLength);

	auto outputIndex = outputOffset;
	for (auto outputIt = output.begin(); outputIt != output.end(); ++outputIt, outputIndex += 1) {
		const auto inputIndex = impl::ChangeSampleRate(sampleRates.Denominator(), sampleRates.Numerator(), outputIndex);
		const auto [firstSampleLoc, secondSampleLoc] = impl::InputIndex2Sample(inputIndex, filter.num_phases());
		const auto firstSampleVal = impl::DotProductSample(input, filter[firstSampleLoc.phaseIndex], firstSampleLoc.inputIndex);
		const auto secondSampleVal = impl::DotProductSample(input, filter[secondSampleLoc.phaseIndex], secondSampleLoc.inputIndex);
		using CommonType = decltype(firstSampleVal);
		*outputIt = (firstSampleVal * CommonType(firstSampleLoc.weight) + secondSampleVal * CommonType(secondSampleLoc.weight))
					/ (CommonType(firstSampleLoc.weight) + CommonType(secondSampleLoc.weight));
	}

	return impl::FindResampleSuspensionPoint(outputIndex, filter.size_original(), filter.num_phases(), sampleRates);
}


/// <summary> Resample the signal using a polyphase FIR low-pass filter. </summary>
/// <param name="input"> The signal to resample. </param>
/// <param name="filter"> The polyphase low-pass filter. </param>
/// <param name="sampleRates"> The { input, output } sample rates. </param>
/// <param name="outputOffset"> Controls the starting point of the output. </param>
/// <remarks>
/// <para> You're free to choose the polyphase filter's parameters, but the more phases
///		the smoother the output. The low-pass cutoff frequency can be calculated with the
///		utility function in this module. </para>
/// <para> Imagine the full (not central) resampled output. The <paramref name="outputOffset"/>
///		gives the fractional index to the first sample of this full output
///		that you want to write into <paramref name="output"/>. </para>
///	</remarks>
template <signal_or_view SignalT,
		  polyphase_or_view PolyphaseTy>
auto Resample(const SignalT& input,
			  const PolyphaseTy& filter,
			  Rational<int64_t> sampleRates,
			  Rational<int64_t> outputOffset,
			  size_t outputLength) {
	using T = typename SignalT::value_type;
	using U = typename PolyphaseTy::container::value_type;
	using R = multiplies_result_t<T, U>;
	constexpr auto Domain = domain_v<SignalT>;

	BasicSignal<R, Domain> out(outputLength, R(0));
	const auto suspensionPoint = Resample(out, input, filter, sampleRates, outputOffset);
	return std::pair{ std::move(out), suspensionPoint };
}


/// <summary> Resample the signal using a polyphase FIR low-pass filter. </summary>
/// <param name="output"> The resampled signal is written here. </param>
/// <param name="input"> The signal to resample. </param>
/// <param name="filter"> The polyphase low-pass filter. </param>
/// <param name="sampleRates"> The { input, output } sample rates. </param>
/// <param name="convMethod"> Full or central convolution. </param>
template <mutable_signal_or_view_r SignalOut,
		  same_domain_as_r<SignalOut> SignalT,
		  polyphase_or_view PolyphaseTy,
		  eConvolutionMethod ConvMethod>
auto Resample(SignalOut&& output,
			  const SignalT& input,
			  const PolyphaseTy& filter,
			  Rational<int64_t> sampleRates,
			  std::integral_constant<eConvolutionMethod, ConvMethod> convMethod) {
	const size_t outputLength = floor(ResampleLength(input.size(), filter.size_original(), filter.num_phases(), sampleRates, convMethod));
	const auto offset = ResampleOffset(input.size(), filter.size_original(), filter.num_phases(), sampleRates, convMethod);
	return Resample(output, input, filter, sampleRates, offset, outputLength);
}


/// <summary> Resample the signal using a polyphase FIR low-pass filter. </summary>
/// <param name="input"> The signal to resample. </param>
/// <param name="filter"> The polyphase low-pass filter. </param>
/// <param name="sampleRates"> The { input, output } sample rates. </param>
/// <param name="convMethod"> Full or central convolution. </param>
template <signal_or_view SignalT,
		  polyphase_or_view PolyphaseTy,
		  eConvolutionMethod ConvMethod>
auto Resample(const SignalT& input,
			  const PolyphaseTy& filter,
			  Rational<int64_t> sampleRates,
			  std::integral_constant<eConvolutionMethod, ConvMethod> convMethod) {
	const size_t outputLength = floor(ResampleLength(input.size(), filter.size_original(), filter.num_phases(), sampleRates, convMethod));
	const auto offset = ResampleOffset(input.size(), filter.size_original(), filter.num_phases(), sampleRates, convMethod);
	return Resample(input, filter, sampleRates, offset, outputLength);
}

} // namespace dspbb