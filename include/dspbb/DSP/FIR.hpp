#include "../Math/Math.hpp"
#include "../Primitives/Signal.hpp"
#include "../Utility/Numbers.hpp"
#include "FFT.hpp"
#include "WindowFunctions.hpp"

#include <numeric>

namespace dspbb {

/// <summary>
/// Creates an arbitrary FIR filter by windowing the impulse response acquired by the IFFT of the frequency response.
/// </summary>
/// <typeparam name="T"> Type of the samples of the returned impulse response. </typeparam>
/// <param name="frequencyResponse"> Amplification of frequencies from 0 to sample rate/2. </param>
/// <param name="numTaps"> Number of coefficients the resulting FIR impulse response.
///		Must be less-equal than twice the number of coefficients in the frequency response. </param>
/// <returns> The coefficients of the impulse response of the FIR filter. </returns>
template <class T>
auto WindowedFirFilter(const Spectrum<T>& frequencyResponse, size_t numTaps = 0) {
	assert(numTaps <= frequencyResponse.Size() * 2);
	Spectrum<T> symmetricResponse(frequencyResponse.Size() * 2);
	std::copy(frequencyResponse.begin(), frequencyResponse.end(), symmetricResponse.begin());
	std::copy(frequencyResponse.rbegin(), frequencyResponse.rend(), symmetricResponse.begin() + frequencyResponse.Size());


	auto impulse = FourierTransform<remove_complex_t<T>>(symmetricResponse);
	auto realImpulse = real(impulse);

	const size_t numCoeffs = realImpulse.Size();
	if (numTaps == 0) {
		numTaps = realImpulse.Size();
	}

	TimeSignal<T> shortImpulse(numTaps);
	std::copy(realImpulse.begin(), realImpulse.begin() + numTaps - numTaps / 2, shortImpulse.begin() + numTaps / 2);
	std::copy(realImpulse.end() - numTaps / 2, realImpulse.end(), shortImpulse.begin());

	return shortImpulse * HammingWindow<T>(numTaps);
}

/// <summary>
/// Creates a low-pass FIR filter by windowing the ideal low-pass filter
/// </summary>
/// <typeparam name="T"> Type of the samples of the returned impulse response. </typeparam>
/// <param name="sampleRate"> The sampling rate of the signal you want to apply the filter to. </param>
/// <param name="cutoffFrequency"> Frequencies below the cutoff are left as is, ones above are silenced. </param>
/// <param name="numTaps"> Number of coefficients of the resulting FIR impulse response. </param>
/// <returns> The coefficients of the impulse response of the FIR LPF. </returns>
template <class T>
auto WindowedLowPass(size_t sampleRate, float cutoffFrequency, size_t numTaps) {
	const T xOffset = T(numTaps - 1) / T(2);
	const T xScale = T(cutoffFrequency) / T(sampleRate) * 2 * pi_v<float>;

	TimeSignal<T> impulse(numTaps);
	for (size_t x = 0; x < numTaps; ++x) {
		T xPrime = xScale * (T(x) - xOffset);
		impulse[x] = std::sin(xPrime) / xPrime;
	}
	if (numTaps % 2 == 1) {
		const size_t centerIndex = (numTaps - 1) / 2;
		impulse[centerIndex] = T(1);
	}

	auto norm = std::accumulate(impulse.begin(), impulse.end(), T(0));
	auto window = HammingWindow<T>(numTaps);
	return impulse * window / norm;
}


/// <summary>
/// Measures how well the <paramref name="filter"/> approximates the desired frequency response.
/// </summary>
/// <param name="filter"> The impulse response of the filter. </param>
/// <param name="desiredResponse"> The desired frequency response of the filter to which the real filter is compared against.
///		Bins range from 0Hz to sampleRate/2, and the number of bins must be at least half of the number of filter taps. </param>
/// <returns> A number between 0 (extremely bad match) and 1 (perfect match). </returns>
template <class T>
auto FilterAccuracy(const TimeSignal<T>& filter, const Spectrum<T>& desiredResponse) {
	if (filter.Size() > 2 * desiredResponse.Size()) {
		throw std::logic_error("You must specify the desired response more accurately with such a large filter.");
	}
	const size_t numBins = desiredResponse.Length() * 2;
	TimeSignal<T> extendedFilter = filter;
	extendedFilter.Resize(numBins, T(0));
	auto actualResponse = abs(FourierTransform(extendedFilter));

	auto magActual = DotProduct(AsConstSpan(actualResponse), AsConstSpan(actualResponse), desiredResponse.Size());
	auto magDesired = DotProduct(AsConstSpan(desiredResponse), AsConstSpan(desiredResponse), desiredResponse.Size());
	auto similarity = DotProduct(AsConstSpan(desiredResponse), AsConstSpan(actualResponse), desiredResponse.Size());

	return similarity / std::max(magActual, magDesired);
}


} // namespace dspbb