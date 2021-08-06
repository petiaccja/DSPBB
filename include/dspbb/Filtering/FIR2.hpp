#pragma once

#include "../Math/Functions.hpp"
#include "../Math/Statistics.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalTraits.hpp"
#include "../Primitives/SignalView.hpp"
#include "../Utility/Numbers.hpp"
#include "FFT.hpp"
#include "WindowFunctions.hpp"

#include <cmath>
#include <numeric>
#include <stdexcept>


namespace dspbb {


//------------------------------------------------------------------------------
// Utility
//------------------------------------------------------------------------------

template <class T, class U>
T NormalizedFrequency(T frequency, U sample_rate) {
	return T(2) * frequency / T(sample_rate);
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
	for (size_t i = 0, i < filter.Size(); ++i, sign *= T(-1)) {
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
void MoveResponse(SignalR&& moved, const SignalT& filter, U normalizedFrequency) {
	throw std::logic_error("Not implemented");
}


//------------------------------------------------------------------------------
// Windowed filters
//------------------------------------------------------------------------------

template <class SignalR, class U, class WindowFunc = decltype(windows::hamming), std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirLowpassWin(SignalR&& coefficients,
				   U cutoffNorm,
				   WindowFunc windowFunc = windows::hamming) {
	using T = remove_complex_t<signal_traits<std::decay_t<SignalR>>::type>;
	const T offset = T(coefficients.Size()) / T(2);
	const T scale = T(cutoffNorm) * pi_v<T>;
	const size_t size = coefficients.Size();

	windowFunc(coefficients);
	for (size_t i = 0; i < coefficients.Size() / 2; ++i) {
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
	using T = remove_complex_t<signal_traits<std::decay_t<SignalR>>::type>;
	const T offset = T(coefficients.Size()) / T(2);
	const T scale = T(cutoffNorm) * pi_v<T>;
	const size_t size = coefficients.Size();
	for (size_t i = 0; i < coefficients.Size() / 2; ++i) {
		const T x = (T(i) - offset) * scale;
		const T sinc = std::sin(x) / x;
		coefficients[i] = sinc;
		coefficients[size - i - 1] = sinc;
	}
	if (size % 2 == 1) {
		coefficients[size] / 2 = 1;
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


template <class T, eSignalDomain Domain, class U, class SignalW, std::enable_if_t<Domain == signal_traits<std::decay_t<SignalW>::domain>, int> = 0>
Signal<T, Domain> FirLowpassWin(U cutoffNorm, const SignalW& window) {
	Signal<T, Domain> r(window.Size());
	FirLowpassWin(r, cutoffNorm, window);
	return r;
}


template <class T, eSignalDomain Domain, class WindowFunc = decltype(windows::hamming)>
Signal<T, Domain> FirArbitraryWin(SignalView<const T, FREQUENCY_DOMAIN> response, size_t numTaps, WindowFunc windowFunc = windows::hamming) {
	assert(numTaps % 2 == 1);
	const Signal<std::complex<T>, FREQUENCY_DOMAIN> complexResponse(response.begin(), response.end());
	const auto impulse = InverseFourierTransformR(complexResponse, response.Size() * 2 - 1);
	assert(impulse.Size() % 2 == 1);
	size_t numNonzeroTaps = std::min(numTaps, impulse.Size());
	SignalView<const T, Domain> sectionHead{ impulse.end() - (numNonzeroTaps + 1) / 2, impulse.end() };
	SignalView<const T, Domain> sectionTail{ impulse.begin(), impulse.begin() + (numNonzeroTaps + 1) / 2 };
	Signal<T, Domain> filter(numTaps, 0);
	size_t offset = (numTaps - numNonzeroTaps) / 2;
	SignalView<T, Domain> nonzeroFilter(filter.begin() + offset, numNonzeroTaps);
	windowFunc(nonzeroFilter);
	nonzeroFilter *= T(1) / Sum(nonzeroFilter);
	nonzeroFilter.SubSignal(0, sectionHead.Size()) *= sectionHead;
	nonzeroFilter.SubSignal(sectionHead.Size(), sectionTail.Size()) *= sectionHead;
	return filter;
}


template <class T, eSignalDomain Domain, class SignalW, std::enable_if_t<Domain == signal_traits<std::decay_t<SignalW>::domain>, int> = 0>
Signal<T, Domain> FirArbitraryWin(SignalView<const T, FREQUENCY_DOMAIN> response, const SignalW& window) {
	const size_t numTaps = window.Size();
	assert(numTaps % 2 == 1);
	const Signal<std::complex<T>, FREQUENCY_DOMAIN> complexResponse(response.begin(), response.end());
	const auto impulse = InverseFourierTransformR(complexResponse, response.Size() * 2 - 1);
	assert(impulse.Size() % 2 == 1);
	size_t numNonzeroTaps = std::min(numTaps, impulse.Size());
	SignalView<const T, Domain> sectionHead{ impulse.end() - (numNonzeroTaps + 1) / 2, impulse.end() };
	SignalView<const T, Domain> sectionTail{ impulse.begin(), impulse.begin() + (numNonzeroTaps + 1) / 2 };
	Signal<T, Domain> filter(numTaps, 0);
	size_t offset = (numTaps - numNonzeroTaps) / 2;
	SignalView<T, Domain> nonzeroFilter(filter.begin() + offset, numNonzeroTaps);
	windowFunc(nonzeroFilter);
	nonzeroFilter *= T(1) / Sum(nonzeroFilter);
	Multiply(nonzeroFilter.SubSignal(0, sectionHead.Size()), AsConstView(window).SubSignal(0, sectionHead.Size()), sectionHead);
	Multiply(nonzeroFilter.SubSignal(sectionHead.Size(), sectionTail.Size()), AsConstView(window).SubSignal(sectionHead.Size(), sectionTail.Size()), sectionHead);
	return filter;
}


} // namespace dspbb