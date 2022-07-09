#pragma once

#include "../../Math/FFT.hpp"
#include "../../Math/Statistics.hpp"
#include "../../Primitives/Signal.hpp"
#include "../../Primitives/SignalView.hpp"
#include "../../Utility/Numbers.hpp"
#include "../Windowing.hpp"


namespace dspbb::fir {

template <class SignalR, class U, class WindowFunc, std::enable_if_t<is_mutable_signal_v<SignalR> && !is_signal_like_v<WindowFunc>, int> = 0>
void KernelWindowedLowpass(SignalR&& coefficients, U cutoffNorm, WindowFunc windowFunc) {
	assert(coefficients.size() % 2 == 1);
	using T = remove_complex_t<typename signal_traits<std::decay_t<SignalR>>::type>;
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


template <class SignalR, class U, class SignalW, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalW>, int> = 0>
void KernelWindowedLowpass(SignalR&& coefficients, U cutoffNorm, const SignalW& window) {
	assert(coefficients.size() % 2 == 1);
	assert(coefficients.size() == window.size());

	using T = remove_complex_t<typename signal_traits<std::decay_t<SignalR>>::type>;
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


template <class SignalR, class ResponseFunc, class WindowFunc, std::enable_if_t<is_mutable_signal_v<SignalR> && std::is_invocable_v<WindowFunc, BasicSignal<float, TIME_DOMAIN>>, int> = 0>
void KernelWindowedArbitrary(SignalR& out, const ResponseFunc& response, WindowFunc windowFunc) {
	assert(out.size() % 2 == 1);
	using R = typename signal_traits<SignalR>::type;
	using ComplexR = std::complex<remove_complex_t<R>>;

	BasicSignal<ComplexR, FREQUENCY_DOMAIN> discreteResponse(out.size() / 2 + 1);
	LinSpace(discreteResponse, R(0), R(1), true);
	std::for_each(discreteResponse.begin(), discreteResponse.end(), [&response](auto& arg) { arg = response(std::real(arg)); });

	const auto impulse = Ifft(discreteResponse, FFT_HALF, out.size() % 2 == 0);
	windowFunc(out);
	AsView(out).subsignal(0, out.size() / 2) *= AsView(impulse).subsignal(impulse.size() / 2 + 1);
	AsView(out).subsignal(out.size() / 2) *= AsView(impulse).subsignal(0, impulse.size() / 2 + 1);
}


template <class SignalR, class ResponseFunc, class SignalW, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalW>, int> = 0>
void KernelWindowedArbitrary(SignalR& out, const ResponseFunc& response, const SignalW& window) {
	assert(out.size() % 2 == 1);
	assert(out.size() == window.size());

	using R = typename signal_traits<SignalR>::type;
	using ComplexR = std::complex<remove_complex_t<R>>;

	BasicSignal<ComplexR, FREQUENCY_DOMAIN> discreteResponse(out.size() / 2 + 1);
	LinSpace(discreteResponse, R(0), R(1), true);
	std::for_each(discreteResponse.begin(), discreteResponse.end(), [&response](auto& arg) { arg = response(std::real(arg)); });

	const auto impulse = Ifft(discreteResponse, FFT_HALF, out.size() % 2 == 0);
	Multiply(AsView(out).subsignal(0, out.size() / 2), AsView(impulse).subsignal(impulse.size() / 2 + 1), AsView(window).subsignal(0, window.size() / 2));
	Multiply(AsView(out).subsignal(out.size() / 2), AsView(impulse).subsignal(0, impulse.size() / 2 + 1), AsView(window).subsignal(window.size() / 2));
}

} // namespace dspbb::fir