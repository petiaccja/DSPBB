#pragma once

#include "../../Math/FFT.hpp"
#include "../../Math/Statistics.hpp"
#include "../../Signal/Signal.hpp"
#include "../../Signal/SignalView.hpp"
#include "../../Signal/Traits.hpp"


namespace dspbb::fir {


/// <summary> Create a low-pass FIR filter using the window method. </summary>
/// <param name="coefficients"> The generated FIR filter. </param>
/// <param name="cutoffNorm"> The normalized cutoff frequency of the low-pass filter. </param>
/// <param name="windowFunc"> The window function factory. </param>
template <mutable_signal_or_view_r SignalR, class U, windows_function_factory WindowFunc>
void KernelWindowedLowpass(SignalR&& coefficients, U cutoffNorm, const WindowFunc& windowFunc) {
	assert(coefficients.size() % 2 == 1);
	using T = remove_complex_t<scalar_type_t<std::decay_t<SignalR>>>;
	const T offset = T(coefficients.size() / 2);
	const T scale = T(cutoffNorm) * pi_v<T>;
	const size_t size = coefficients.size();

	windowFunc(coefficients);
	for (size_t i = 0; i < size / 2; ++i) {
		const T x = (T(i) - offset) * scale;
		const T sinc = std::sin(x) / x;
		coefficients[i] *= sinc;
		coefficients[size - i - 1] *= sinc;
	}
	coefficients *= T(1) / T(Sum(coefficients));
}


/// <summary> Create a low-pass FIR filter using the window method. </summary>
///	<param name="coefficients"> The generated FIR filter. </param>
/// <param name="cutoffNorm"> The normalized cutoff frequency of the low-pass filter. </param>
/// <param name="window"> The coefficients of the window function. </param>
template <mutable_signal_or_view_r SignalR, class U, same_domain_as_r<SignalR> SignalW>
void KernelWindowedLowpass(SignalR&& coefficients, U cutoffNorm, const SignalW& window) {
	assert(!IsAliasing(coefficients, window));
	assert(coefficients.size() % 2 == 1);
	assert(coefficients.size() == window.size());

	using T = remove_complex_t<scalar_type_t<std::decay_t<SignalR>>>;
	const T offset = T(coefficients.size() / 2);
	const T scale = T(cutoffNorm) * pi_v<T>;
	const size_t size = coefficients.size();
	for (size_t i = 0; i < size / 2; ++i) {
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


/// <summary> Create a low-pass FIR filter using the window method. </summary>
/// <param name="coefficients"> The generated FIR filter. </param>
/// <param name="response"> The continuous response of the filter. </param>
/// <param name="windowFunc"> The window function factory. </param>
template <mutable_signal_or_view_r SignalR, class ResponseFunc, windows_function_factory WindowFunc>
void KernelWindowedArbitrary(SignalR& coefficients, const ResponseFunc& response, const WindowFunc& windowFunc) {
	assert(coefficients.size() % 2 == 1);
	using R = scalar_type_t<SignalR>;
	using ComplexR = std::complex<remove_complex_t<R>>;

	BasicSignal<ComplexR, FREQUENCY_DOMAIN> discreteResponse(coefficients.size() / 2 + 1);
	LinSpace(discreteResponse, R(0), R(1), true);
	std::for_each(discreteResponse.begin(), discreteResponse.end(), [&response](auto& arg) { arg = response(std::real(arg)); });

	const auto impulse = Ifft(discreteResponse, FFT_HALF, coefficients.size() % 2 == 0);
	windowFunc(coefficients);
	AsView(coefficients).subsignal(0, coefficients.size() / 2) *= AsView(impulse).subsignal(impulse.size() / 2 + 1);
	AsView(coefficients).subsignal(coefficients.size() / 2) *= AsView(impulse).subsignal(0, impulse.size() / 2 + 1);
}


/// <summary> Create a low-pass FIR filter using the window method. </summary>
/// <param name="coefficients"> The generated FIR filter. </param>
/// <param name="response"> The continuous response of the filter. </param>
/// <param name="window"> The coefficients of the window function. </param>
template <mutable_signal_or_view_r SignalR, class ResponseFunc, same_domain_as_r<SignalR> SignalW>
void KernelWindowedArbitrary(SignalR& coefficients, const ResponseFunc& response, const SignalW& window) {
	assert(coefficients.size() % 2 == 1);
	assert(coefficients.size() == window.size());

	using R = scalar_type_t<SignalR>;
	using ComplexR = std::complex<remove_complex_t<R>>;

	BasicSignal<ComplexR, FREQUENCY_DOMAIN> discreteResponse(coefficients.size() / 2 + 1);
	LinSpace(discreteResponse, R(0), R(1), true);
	std::for_each(discreteResponse.begin(), discreteResponse.end(), [&response](auto& arg) { arg = response(std::real(arg)); });

	const auto impulse = Ifft(discreteResponse, FFT_HALF, coefficients.size() % 2 == 0);
	Multiply(AsView(coefficients).subsignal(0, coefficients.size() / 2), AsView(impulse).subsignal(impulse.size() / 2 + 1), AsView(window).subsignal(0, window.size() / 2));
	Multiply(AsView(coefficients).subsignal(coefficients.size() / 2), AsView(impulse).subsignal(0, impulse.size() / 2 + 1), AsView(window).subsignal(window.size() / 2));
}

} // namespace dspbb::fir