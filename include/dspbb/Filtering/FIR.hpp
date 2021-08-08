#pragma once

#include "../Math/Functions.hpp"
#include "../Math/Statistics.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "../Utility/Numbers.hpp"
#include "FFT.hpp"
#include "WindowFunctions.hpp"


namespace dspbb {


//------------------------------------------------------------------------------
// Utility
//------------------------------------------------------------------------------

template <class T, class U>
T NormalizedFrequency(T frequency, U sampleRate) {
	return T(2) * frequency / T(sampleRate);
}


//------------------------------------------------------------------------------
// Band transforms
//------------------------------------------------------------------------------

template <class SignalR, class SignalT, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
void MirrorResponse(SignalR&& mirrored, const SignalT& filter) {
	assert(mirrored.Size() == filter.Size());
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	T sign = T(1);
	for (size_t i = 0; i < filter.Size(); ++i, sign *= T(-1)) {
		mirrored[i] = R(sign * filter[i]);
	}
}

template <class SignalR, class SignalT, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
void ComplementaryResponse(SignalR&& complementary, const SignalT& filter) {
	assert(filter.Size() % 2 == 1);
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	Multiply(complementary, filter, T(-1));
	complementary[complementary.Size() / 2] += R(1);
}

template <class SignalR, class SignalT, class U, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
void ShiftResponse(SignalR&& moved, const SignalT& filter, U normalizedFrequency) {
	U offset = static_cast<U>(filter.Size() / 2);
	U scale = pi_v<U> * normalizedFrequency;
	const size_t size = filter.Size();
	for (size_t i = 0; i < size / 2; ++i) {
		const U x = (U(i) - offset) * scale;
		const U c = std::cos(x);
		moved[i] = c * filter[i];
		moved[size - i - 1] = c * filter[size - i - 1];
	}
	moved *= typename signal_traits<SignalT>::type(2);
}


//------------------------------------------------------------------------------
// Windowed filters
//------------------------------------------------------------------------------

//--------------------------------------
// Lowpass
//--------------------------------------

template <class SignalR, class U, class WindowFunc = decltype(windows::hamming), std::enable_if_t<is_mutable_signal_v<SignalR> && !is_signal_like_v<WindowFunc>, int> = 0>
void FirLowpassWin(SignalR&& coefficients,
				   U cutoffNorm,
				   WindowFunc windowFunc = windows::hamming) {
	assert(coefficients.Size() % 2 == 1);
	using T = remove_complex_t<signal_traits<std::decay_t<SignalR>>::type>;
	const T offset = T(coefficients.Size() / 2);
	const T scale = T(cutoffNorm) * pi_v<T>;
	const size_t size = coefficients.Size();

	windowFunc(coefficients);
	for (size_t i = 0; i < size / 2; ++i) {
		const T x = (T(i) - offset) * scale;
		const T sinc = std::sin(x) / x;
		coefficients[i] *= sinc;
		coefficients[size - i - 1] *= sinc;
	}
	coefficients *= T(1) / T(Sum(coefficients));
}


template <class SignalR, class U, class SignalW, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalW>, int> = 0>
void FirLowpassWin(SignalR&& coefficients,
				   U cutoffNorm,
				   const SignalW& window) {
	assert(coefficients.Size() % 2 == 1);
	using T = remove_complex_t<signal_traits<std::decay_t<SignalR>>::type>;
	const T offset = T(coefficients.Size() / 2);
	const T scale = T(cutoffNorm) * pi_v<T>;
	const size_t size = coefficients.Size();
	for (size_t i = 0; i < coefficients.Size() / 2; ++i) {
		const T x = (T(i) - offset) * scale;
		const T sinc = std::sin(x) / x;
		coefficients[i] = sinc;
		coefficients[size - i - 1] = sinc;
	}
	if (size % 2 == 1) {
		coefficients[size / 2] = 1;
	}
	coefficients *= window;
	coefficients *= T(1) / T(Sum(coefficients));
}


template <class T, eSignalDomain Domain, class U, class WindowFunc = decltype(windows::hamming)>
Signal<T, Domain> FirLowpassWin(U cutoffNorm,
								size_t numTaps,
								WindowFunc windowFunc = windows::hamming) {
	Signal<T, Domain> r(numTaps);
	FirLowpassWin(r, cutoffNorm, windowFunc);
	return r;
}


template <class T, eSignalDomain Domain, class U, class SignalW, std::enable_if_t<Domain == signal_traits<std::decay_t<SignalW>>::domain, int> = 0>
Signal<T, Domain> FirLowpassWin(U cutoffNorm,
								const SignalW& window) {
	Signal<T, Domain> r(window.Size());
	FirLowpassWin(r, cutoffNorm, window);
	return r;
}


//--------------------------------------
// Highpass
//--------------------------------------

template <class SignalR, class U, class WindowFunc = decltype(windows::hamming), std::enable_if_t<is_mutable_signal_v<SignalR> && !is_signal_like_v<WindowFunc>, int> = 0>
void FirHighpassWin(SignalR&& coefficients,
					U cutoffNorm,
					WindowFunc windowFunc = windows::hamming) {
	FirLowpassWin(coefficients, cutoffNorm, windowFunc);
	ComplementaryResponse(coefficients, coefficients);
}
template <class SignalR, class U, class SignalW, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalW>, int> = 0>
void FirHighpassWin(SignalR&& coefficients,
					U cutoffNorm,
					const SignalW& window) {
	FirLowpassWin(coefficients, cutoffNorm, window);
	ComplementaryResponse(coefficients, coefficients);
}
template <class T, eSignalDomain Domain, class U, class WindowFunc = decltype(windows::hamming)>
Signal<T, Domain> FirHighpassWin(U cutoffNorm,
								 size_t numTaps,
								 WindowFunc windowFunc = windows::hamming) {
	Signal<T, Domain> r(numTaps);
	FirHighpassWin(r, cutoffNorm, windowFunc);
	return r;
}
template <class T, eSignalDomain Domain, class U, class SignalW, std::enable_if_t<Domain == signal_traits<std::decay_t<SignalW>>::domain, int> = 0>
Signal<T, Domain> FirHighpassWin(U cutoffNorm,
								 const SignalW& window) {
	Signal<T, Domain> r(window.Size());
	FirHighpassWin(r, cutoffNorm, window);
	return r;
}


//--------------------------------------
// Bandpass
//--------------------------------------

template <class SignalR, class U, class WindowFunc = decltype(windows::hamming), std::enable_if_t<is_mutable_signal_v<SignalR> && !is_signal_like_v<WindowFunc>, int> = 0>
void FirBandpassWin(SignalR&& coefficients,
					U bandLow,
					U bandHigh,
					WindowFunc windowFunc = windows::hamming) {
	U bandWidth = bandHigh - bandLow;
	U bandCenter = (bandHigh + bandLow) / 2;
	FirLowpassWin(coefficients, bandWidth / U(2), windowFunc);
	ShiftResponse(coefficients, coefficients, bandCenter);
}
template <class SignalR, class U, class SignalW, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalW>, int> = 0>
void FirBandpassWin(SignalR&& coefficients,
					U bandLow,
					U bandHigh,
					const SignalW& window) {
	U bandWidth = bandHigh - bandLow;
	U bandCenter = (bandHigh + bandLow) / 2;
	FirLowpassWin(coefficients, bandWidth / U(2), window);
	ShiftResponse(coefficients, coefficients, bandCenter);
}
template <class T, eSignalDomain Domain, class U, class WindowFunc = decltype(windows::hamming)>
Signal<T, Domain> FirBandpassWin(U bandLow,
								 U bandHigh,
								 size_t numTaps,
								 WindowFunc windowFunc = windows::hamming) {
	Signal<T, Domain> r(numTaps);
	FirBandpassWin(r, bandLow, bandHigh, windowFunc);
	return r;
}
template <class T, eSignalDomain Domain, class U, class SignalW, std::enable_if_t<Domain == signal_traits<std::decay_t<SignalW>>::domain, int> = 0>
Signal<T, Domain> FirBandpassWin(U bandLow,
								 U bandHigh,
								 const SignalW& window) {
	Signal<T, Domain> r(window.Size());
	FirBandpassWin(r, bandLow, bandHigh, window);
	return r;
}


//--------------------------------------
// Arbitrary response
//--------------------------------------

template <class T, eSignalDomain Domain, class WindowFunc = decltype(windows::hamming)>
Signal<T, Domain> FirArbitraryWin(SignalView<const T, FREQUENCY_DOMAIN> response, size_t numTaps, WindowFunc windowFunc = windows::hamming) {
	assert(numTaps % 2 == 1);
	const Signal<std::complex<T>, FREQUENCY_DOMAIN> complexResponse(response.begin(), response.end());
	const auto impulse = InverseFourierTransformR(complexResponse, response.Size() * 2 - 1);
	assert(impulse.Size() % 2 == 1);
	size_t numNonzeroTaps = std::min(numTaps, impulse.Size());
	SignalView<const T, Domain> sectionHead{ impulse.end() - numNonzeroTaps / 2, impulse.end() };
	SignalView<const T, Domain> sectionTail{ impulse.begin(), impulse.begin() + (numNonzeroTaps + 1) / 2 };
	Signal<T, Domain> filter(numTaps, T(0));
	size_t offset = (numTaps - numNonzeroTaps) / 2;
	SignalView<T, Domain> nonzeroFilter(filter.begin() + offset, numNonzeroTaps);
	windowFunc(nonzeroFilter);
	nonzeroFilter.SubSignal(0, sectionHead.Size()) *= sectionHead;
	nonzeroFilter.SubSignal(sectionHead.Size(), sectionTail.Size()) *= sectionTail;
	return filter;
}


template <class T, eSignalDomain Domain, class SignalW, std::enable_if_t<Domain == signal_traits<std::decay_t<SignalW>>::domain, int> = 0>
Signal<T, Domain> FirArbitraryWin(SignalView<const T, FREQUENCY_DOMAIN> response, const SignalW& window) {
	const size_t numTaps = window.Size();
	assert(numTaps % 2 == 1);
	const Signal<std::complex<T>, FREQUENCY_DOMAIN> complexResponse(response.begin(), response.end());
	const auto impulse = InverseFourierTransformR(complexResponse, response.Size() * 2 - 1);
	assert(impulse.Size() % 2 == 1);
	size_t numNonzeroTaps = std::min(numTaps, impulse.Size());
	SignalView<const T, Domain> sectionHead{ impulse.end() - numNonzeroTaps / 2, impulse.end() };
	SignalView<const T, Domain> sectionTail{ impulse.begin(), impulse.begin() + (numNonzeroTaps + 1) / 2 };
	Signal<T, Domain> filter(numTaps, T(0));
	size_t offset = (numTaps - numNonzeroTaps) / 2;
	SignalView<T, Domain> nonzeroFilter(filter.begin() + offset, numNonzeroTaps);
	Multiply(nonzeroFilter.SubSignal(0, sectionHead.Size()), AsConstView(window).SubSignal(0, sectionHead.Size()), sectionHead);
	Multiply(nonzeroFilter.SubSignal(sectionHead.Size(), sectionTail.Size()), AsConstView(window).SubSignal(sectionHead.Size(), sectionTail.Size()), sectionTail);
	return filter;
}


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