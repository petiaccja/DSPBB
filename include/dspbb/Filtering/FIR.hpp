#pragma once

#include "../Math/Functions.hpp"
#include "../Math/Statistics.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "../Utility/Numbers.hpp"
#include "FFT.hpp"
#include "WindowFunctions.hpp"

#include <numeric>

namespace dspbb {

template <class T>
struct FirGeneralQuality {
	T cosineSimilarity;
	T maxDifference;
	T maxRelDifference;
	T meanDifference;
	T meanRelDifference;
	T meanSquareDifference;
	T meanSquareRelDifference;
};

/// <summary>
/// Creates an arbitrary FIR filter by windowing the impulse response acquired by the IFFT of the frequency response.
/// </summary>
/// <typeparam name="T"> Type of the samples of the returned impulse response. </typeparam>
/// <param name="frequencyResponse"> Amplification of frequencies from 0 to sample rate/2. </param>
/// <param name="window"> The window function to apply to the impulse response. Also defines number of taps,
///		which must be less-equal than twice the number of coefficients in the frequency response. </param>
/// <returns> The coefficients of the impulse response of the FIR filter. </returns>
template <class T>
TimeSignal<T> FirGeneralWindowed(const Spectrum<T>& frequencyResponse, TimeSignalView<const T> window) {
	const size_t numTaps = window.Size();
	assert(numTaps != 0);

	auto impulse = InverseFourierTransformR(Spectrum<std::complex<T>>(frequencyResponse.begin(), frequencyResponse.end()), frequencyResponse.Size() * 2 - 1);
	assert(numTaps <= impulse.Size());
	auto realImpulse = Real(impulse);

	TimeSignal<T> shortImpulse(numTaps);
	std::copy(realImpulse.begin(), realImpulse.begin() + numTaps - numTaps / 2, shortImpulse.begin() + numTaps / 2);
	std::copy(realImpulse.end() - numTaps / 2, realImpulse.end(), shortImpulse.begin());

	return shortImpulse * window;
}


/// <summary>
/// Creates an arbitrary FIR filter by windowing the impulse response acquired by the IFFT of the frequency response.
/// </summary>
/// <typeparam name="T"> Type of the samples of the returned impulse response. </typeparam>
/// <param name="frequencyResponse"> Amplification of frequencies from 0 to sample rate/2. </param>
/// <param name="window"> The window function to apply to the impulse response. Also defines number of taps,
///		which must be less-equal than twice the number of coefficients in the frequency response. </param>
/// <returns> The coefficients of the impulse response of the FIR filter. </returns>
template <class T>
TimeSignal<T> FirGeneralWindowed(const Spectrum<T>& frequencyResponse, const TimeSignal<T>& window) {
	return FirGeneralWindowed(frequencyResponse, AsConstView(window));
}


/// <summary>
/// Creates an arbitrary FIR filter by windowing the impulse response acquired by the IFFT of the frequency response.
/// </summary>
/// <typeparam name="T"> Type of the samples of the returned impulse response. </typeparam>
/// <param name="frequencyResponse"> Amplification of frequencies from 0 to sample rate/2. </param>
/// <param name="numTaps"> Number of coefficients the resulting FIR impulse response.
///		Must be less-equal than twice the number of coefficients in the frequency response. </param>
/// <returns> The coefficients of the impulse response of the FIR filter. </returns>
template <class T, class WindowFunc = decltype(windows::hamming)>
TimeSignal<T> FirGeneralWindowed(const Spectrum<T>& frequencyResponse,
								 size_t numTaps,
								 WindowFunc windowFunc = windows::hamming) {
	return FirGeneralWindowed(frequencyResponse, windowFunc.operator()<T, TIME_DOMAIN>(numTaps));
}


/// <summary>
/// Creates a low-pass FIR filter by windowing the ideal low-pass filter.
/// </summary>
/// <typeparam name="T"> Type of the samples of the returned impulse response. </typeparam>
/// <param name="sampleRate"> The sampling rate of the signal you want to apply the filter to. </param>
/// <param name="cutoffFrequency"> Frequencies below the cutoff are left as is, ones above are silenced. </param>
/// <param name="window"> The window to apply to the ideal filter to reduce ringing. Also defines the number of taps. </param>
/// <returns></returns>
template <class T>
TimeSignal<T> FirLowPassWindowed(float cutoffFrequency, size_t sampleRate, TimeSignalView<const T> window) {
	const auto numTaps = window.Size();
	const T xOffset = T(numTaps - 1) / T(2);
	const T xScale = T(cutoffFrequency) / T(sampleRate) * T(2) * pi_v<T>;

	TimeSignal<T> impulse(numTaps);
	for (size_t x = 0; x < numTaps; ++x) {
		T xPrime = xScale * (T(x) - xOffset);
		impulse[x] = std::sin(xPrime) / xPrime;
	}
	if (numTaps % 2 == 1) {
		const size_t centerIndex = (numTaps - 1) / 2;
		impulse[centerIndex] = T(1);
	}

	impulse *= window;
	impulse *= T(1) / Sum(impulse);
	return impulse;
}


/// <summary>
/// Creates a low-pass FIR filter by windowing the ideal low-pass filter.
/// </summary>
/// <typeparam name="T"> Type of the samples of the returned impulse response. </typeparam>
/// <param name="sampleRate"> The sampling rate of the signal you want to apply the filter to. </param>
/// <param name="cutoffFrequency"> Frequencies below the cutoff are left as is, ones above are silenced. </param>
/// <param name="window"> The window to apply to the ideal filter to reduce ringing. Also defines the number of taps. </param>
/// <returns></returns>
template <class T>
TimeSignal<T> FirLowPassWindowed(float cutoffFrequency, size_t sampleRate, const TimeSignal<T>& window) {
	return FirLowPassWindowed(cutoffFrequency, sampleRate, AsConstView(window));
}


/// <summary>
/// Creates a low-pass FIR filter by windowing the ideal low-pass filter.
/// </summary>
/// <typeparam name="T"> Type of the samples of the returned impulse response. </typeparam>
/// <param name="sampleRate"> The sampling rate of the signal you want to apply the filter to. </param>
/// <param name="cutoffFrequency"> Frequencies below the cutoff are left as is, ones above are silenced. </param>
/// <param name="numTaps"> Number of coefficients of the resulting FIR impulse response. </param>
/// <returns> The coefficients of the impulse response of the FIR LPF. </returns>
template <class T, class WindowFunc = decltype(windows::hamming)>
TimeSignal<T> FirLowPassWindowed(float cutoffFrequency,
								 size_t sampleRate,
								 size_t numTaps,
								 WindowFunc windowFunc = windows::hamming) {
	return FirLowPassWindowed(cutoffFrequency, sampleRate, windowFunc.operator()<T, TIME_DOMAIN>(numTaps));
}


/// <summary>
/// Measures how well the <paramref name="filter"/> approximates the desired frequency response.
/// </summary>
/// <param name="filter"> The impulse response of the filter. </param>
/// <param name="desiredResponse"> The desired frequency response of the filter to which the real filter is compared against.
///		Bins range from 0Hz to sampleRate/2, and the number of bins must be at least half of the number of filter taps. </param>
/// <returns> A number between 0 (extremely bad match) and 1 (perfect match). </returns>
template <class T>
FirGeneralQuality<T> FirQuality(const TimeSignal<T>& filter, const Spectrum<T>& desiredResponse) {
	const size_t numBins = 2 * desiredResponse.Length() - 1 - desiredResponse.Length() % 2;
	assert(numBins >= filter.Size());
	TimeSignal<T> extendedFilter = filter;
	extendedFilter.Resize(numBins, T(0));
	const auto actualResponse = Abs(FourierTransform(extendedFilter, false));
	assert(actualResponse.Size() == desiredResponse.Size());
	const auto difference = actualResponse - desiredResponse;
	const auto relDifference = difference / desiredResponse;

	auto magActual = DotProduct(AsConstView(actualResponse), AsConstView(actualResponse), desiredResponse.Size());
	auto magDesired = DotProduct(AsConstView(desiredResponse), AsConstView(desiredResponse), desiredResponse.Size());
	auto similarity = DotProduct(AsConstView(desiredResponse), AsConstView(actualResponse), desiredResponse.Size());

	FirGeneralQuality<T> quality;
	quality.cosineSimilarity = similarity / std::max(magActual, magDesired);
	quality.maxDifference = Max(Abs(difference));
	quality.maxRelDifference = Max(Abs(relDifference));
	quality.meanDifference = Mean(Abs(difference));
	quality.meanRelDifference = Mean(Abs(relDifference));
	quality.meanSquareDifference = RootMeanSquare(difference);
	quality.meanSquareRelDifference = RootMeanSquare(relDifference);

	return quality;
}


} // namespace dspbb