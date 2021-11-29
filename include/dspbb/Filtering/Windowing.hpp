#pragma once

#include "../Generators/Spaces.hpp"
#include "../Math/FFT.hpp"
#include "../Math/Statistics.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "../Utility/Numbers.hpp"
#include "../Utility/TypeTraits.hpp"

#include <cmath>


namespace dspbb {

//------------------------------------------------------------------------------
// Assess properties of windows.
//------------------------------------------------------------------------------
template <class T, eSignalDomain Domain>
T CoherentGain(BasicSignalView<T, Domain> window) {
	return Sum(window) / remove_complex_t<T>(window.Size());
}

template <class T, eSignalDomain Domain>
T CoherentGain(const BasicSignal<T, Domain>& window) {
	return CoherentGain(AsConstView(window));
}

template <class T, eSignalDomain Domain>
T EnergyGain(BasicSignalView<T, Domain> window) {
	return SumSquare(window) / T(window.Size());
}

template <class T, eSignalDomain Domain>
T EnergyGain(const BasicSignal<T, Domain>& window) {
	return EnergyGain(AsConstView(window));
}


//------------------------------------------------------------------------------
// List of window functions.
//------------------------------------------------------------------------------
template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void HammingWindow(SignalR&& out) {
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using U = remove_complex_t<R>;

	LinSpace(out, U(0), U(2) * pi_v<U>, true);
	Cos(out, out);
	out *= U(-0.46);
	out += U(0.54);
}

template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FlatTopWindow(SignalR&& out) {
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using U = remove_complex_t<R>;

	U c0 = U(0.21557895);
	U c1 = U(-0.41663158);
	U c2 = U(0.277263158);
	U c3 = U(-0.083578947);
	U c4 = U(0.006947368);

	const U N = U(out.Length());
	U preSize1 = U(2) * pi_v<U> / (N - U(1));
	U preSize2 = U(4) * pi_v<U> / (N - U(1));
	U preSize3 = U(6) * pi_v<U> / (N - U(1));
	U preSize4 = U(8) * pi_v<U> / (N - U(1));

	LinSpace(out, U(0), U(out.Size() - 1), true);
	std::for_each(out.begin(), out.end(), [&](R& k) {
		const U kreal = std::real(k);
		k = c0
			+ c1 * std::cos(preSize1 * kreal)
			+ c2 * std::cos(preSize2 * kreal)
			+ c3 * std::cos(preSize3 * kreal)
			+ c4 * std::cos(preSize4 * kreal);
	});
}

template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void RectangularWindow(SignalR&& out) {
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using U = remove_complex_t<R>;
	std::fill(out.begin(), out.end(), R(U(1.0)));
}

template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void TriangularWindow(SignalR&& out) {
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using U = remove_complex_t<R>;
	LinSpace(out, U(0), U(2), true);
	out -= U(1);
	Abs(out, out);
	out *= U(-1);
	out += U(1);
}

template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void BlackmanWindow(SignalR&& out) {
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using U = remove_complex_t<R>;
	LinSpace(out, U(0), U(2) * pi_v<U>, true);
	std::for_each(out.begin(), out.end(), [&](R& k) {
		const U kreal = std::real(k);
		k = U(0.42) - U(0.5) * std::cos(kreal) + U(0.08) * std::cos(2 * kreal);
	});
}

template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void BlackmanHarrisWindow(SignalR&& out) {
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using U = remove_complex_t<R>;
	LinSpace(out, U(0), U(2) * pi_v<U>, true);
	std::for_each(out.begin(), out.end(), [&](R& k) {
		const U kreal = std::real(k);
		k = U(0.35875) - U(0.48829) * std::cos(kreal) + U(0.14128) * std::cos(2 * kreal) + U(-0.01168) * std::cos(3 * kreal);
	});
}

template <class SignalR, class V, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void GaussianWindow(SignalR&& out, V sigma = 1.f) {
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using U = remove_complex_t<R>;
	const auto N = U(out.Size());
	const auto M = (N - U(1)) / U(2);
	LinSpace(out, -M, M, true);
	out *= U(1) / (U(sigma) * M);
	Multiply(out, out, out);
	out *= U(-0.5);
	Exp(out, out);
}

template <class SignalR, class V, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void KaiserWindow(SignalR&& out, V alpha) {
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using U = remove_complex_t<R>;
	LinSpace(out, -U(1), U(1), true);
	std::for_each(out.begin(), out.end(), [&](R& k) {
		const U kreal = std::real(k);
		const U piAlpha = pi_v<U> * U(alpha);
		const U arg = std::sqrt(std::max(U(0), U(1) - kreal * kreal));
		k = U(std::cyl_bessel_i(U(0), piAlpha * arg)) / U(std::cyl_bessel_i(U(0), U(piAlpha)));
	});
}

template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void LanczosWindow(SignalR&& out) {
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using U = remove_complex_t<R>;
	LinSpace(out, -pi_v<U>, pi_v<U>, true);
	std::for_each(out.begin(), out.end(), [&](R& k) {
		const U kreal = std::real(k);
		k = kreal != U(0) ? std::sin(kreal) / kreal : U(1);
	});
}

namespace impl {
	template <class T>
	auto ChebyshevPoly(size_t n, const T& x) {
		const T tn = static_cast<T>(n);
		const T sign = n % 2 == 0 ? T(1) : T(-1);
		if (x < -T(1)) {
			return sign * std::cosh(tn * std::acosh(-x));
		}
		if (x <= T(1)) {
			return std::cos(tn * std::acos(x));
		}
		return std::cosh(tn * std::acosh(x));
	};
} // namespace impl

template <class SignalR, class V, std::enable_if_t<is_mutable_signal_v<SignalR> && !is_complex_v<typename std::decay_t<SignalR>::value_type>, int> = 0>
void DolphChebyshevWindow(SignalR&& out, V attenuation) {
	using T = typename std::decay_t<SignalR>::value_type;

	const size_t M = out.Size() - 1;
	const T beta = std::cosh(T(1) / M * std::acosh(T(1) / T(attenuation)));
	Spectrum<std::complex<T>> spectrum(out.Size() / 2 + 1);
	LinSpace(spectrum, T(0), pi_v<T> * (T(spectrum.Size() - 1) / T(out.Size())), true);
	std::for_each(spectrum.begin(), spectrum.end(), [M, beta](auto& k) {
		const auto i = std::complex<T>(0, 1);
		const auto phase = std::exp(i * k * T(M % 2));
		const auto amplitude = impl::ChebyshevPoly(M, beta * std::cos(std::real(k)));
		k = phase * amplitude;
	});

	Ifft(out, spectrum);
	FftShift(out, out);
	const T normalization = kernels::MapReduceVectorized(
		out.Data(), out.Size(), T(0),
		[](const auto& acc, const auto& v) { return kernels::math_functions::max(acc, v); },
		[](const auto& v) { return kernels::math_functions::abs(v); });
	out *= T(1) / normalization;
}

template <class SignalR, class V, std::enable_if_t<is_mutable_signal_v<SignalR> && is_complex_v<typename std::decay_t<SignalR>::value_type>, int> = 0>
void DolphChebyshevWindow(SignalR&& out, V attenuation) {
	using R = typename std::decay_t<SignalR>::value_type;
	using T = remove_complex_t<R>;
	constexpr auto domain = signal_traits<std::decay_t<SignalR>>::domain;

	BasicSignal<T, domain> outReal(out.Size());
	DolphChebyshevWindow(outReal, attenuation);
	std::transform(outReal.begin(), outReal.end(), out.begin(), [](auto& v) { return R{ v, T(0) }; });
}


template <class T, eSignalDomain Domain = eSignalDomain::TIME>
BasicSignal<T, Domain> HammingWindow(size_t length) {
	BasicSignal<T, Domain> window(length);
	HammingWindow(AsView(window));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
BasicSignal<T, Domain> FlatTopWindow(size_t length) {
	BasicSignal<T, Domain> window(length);
	FlatTopWindow(AsView(window));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
BasicSignal<T, Domain> RectangularWindow(size_t length) {
	BasicSignal<T, Domain> window(length, T(1.0));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
BasicSignal<T, Domain> TriangularWindow(size_t length) {
	BasicSignal<T, Domain> window(length);
	TriangularWindow(AsView(window));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
BasicSignal<T, Domain> BlackmanWindow(size_t length) {
	BasicSignal<T, Domain> window(length);
	BlackmanWindow(AsView(window));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
BasicSignal<T, Domain> BlackmanHarrisWindow(size_t length) {
	BasicSignal<T, Domain> window(length);
	BlackmanHarrisWindow(AsView(window));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
BasicSignal<T, Domain> GaussianWindow(size_t length, T sigma = T(1)) {
	BasicSignal<T, Domain> window(length);
	GaussianWindow(AsView(window), sigma);
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
BasicSignal<T, Domain> KaiserWindow(size_t length, T alpha = T(1)) {
	BasicSignal<T, Domain> window(length);
	KaiserWindow(AsView(window), alpha);
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
BasicSignal<T, Domain> LanczosWindow(size_t length) {
	BasicSignal<T, Domain> window(length);
	LanczosWindow(AsView(window));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
BasicSignal<T, Domain> DolphChebyshevWindow(size_t length, T attenuation) {
	BasicSignal<T, Domain> window(length);
	DolphChebyshevWindow(AsView(window), attenuation);
	return window;
}

//------------------------------------------------------------------------------
// Helper for when you have to pass a window function as an argument.
//------------------------------------------------------------------------------
namespace windows {
	struct Hamming {
		template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
		auto operator()(SignalR&& out) const {
			return HammingWindow(out);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return HammingWindow<T, Domain>(length);
		}
	} const hamming;

	struct Rectangular {
		template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
		auto operator()(SignalR&& out) const {
			return RectangularWindow(out);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return RectangularWindow<T, Domain>(length);
		}
	} const rectangular;

	struct Flattop {
		template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
		auto operator()(SignalR&& out) const {
			return FlatTopWindow(out);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return FlatTopWindow<T, Domain>(length);
		}
	} const flattop;

	struct Triangular {
		template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
		auto operator()(SignalR&& out) const {
			return TriangularWindow(out);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return TriangularWindow<T, Domain>(length);
		}
	} const triangular;

	struct Blackman {
		template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
		auto operator()(SignalR&& out) const {
			return BlackmanWindow(out);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return BlackmanWindow<T, Domain>(length);
		}
	} const blackman;

	struct BlackmanHarris {
		template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
		auto operator()(SignalR&& out) const {
			return BlackmanHarrisWindow(out);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return BlackmanHarrisWindow<T, Domain>(length);
		}
	} const blackmanHarris;

	struct Gaussian {
		template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
		auto operator()(SignalR&& out) const {
			return GaussianWindow(out, m_sigma);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return GaussianWindow<T, Domain>(length, T(m_sigma));
		}
		template <class T>
		Gaussian sigma(T sigma) const {
			auto copy = *this;
			copy.m_sigma = double(sigma);
			return copy;
		}
		double m_sigma = 1;
	} const gaussian;

	struct Kaiser {
		template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
		auto operator()(SignalR&& out) const {
			return KaiserWindow(out, m_alpha);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return KaiserWindow<T, Domain>(length, T(m_alpha));
		}
		template <class T>
		Kaiser alpha(T alpha) const {
			auto copy = *this;
			copy.m_alpha = double(alpha);
			return copy;
		}
		double m_alpha = 1;
	} const kaiser;

	struct Lanczos {
		template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
		auto operator()(SignalR&& out) const {
			return LanczosWindow(out);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return LanczosWindow<T, Domain>(length);
		}
	} const lanczos;

	struct DolphChebyshev {
		template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
		auto operator()(SignalR&& out) const {
			return DolphChebyshevWindow(out, m_attenuation);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return DolphChebyshevWindow<T, Domain>(length, T(m_attenuation));
		}
		template <class T>
		DolphChebyshev attenuation(T atten) const {
			auto copy = *this;
			copy.m_attenuation = double(atten);
			return copy;
		}
		double m_attenuation = 1;
	} const dolphChebyshev;

} // namespace windows


} // namespace dspbb