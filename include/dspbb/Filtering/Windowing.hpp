#pragma once

#include "../Generators/Spaces.hpp"
#include "../Math/FFT.hpp"
#include "../Math/Statistics.hpp"
#include "../Signal/Signal.hpp"
#include "../Signal/SignalView.hpp"
#include "../Utility/Numbers.hpp"
#include "../Utility/TypeTraits.hpp"

#include <cmath>


namespace dspbb {

//------------------------------------------------------------------------------
// Assess properties of windows.
//------------------------------------------------------------------------------

template <signal_or_view SignalTy>
auto CoherentGain(const SignalTy& window) {
	using Scalar = remove_complex_t<typename SignalTy::value_type>;
	return Sum(window) / Scalar(window.size());
}


template <signal_or_view SignalTy>
auto EnergyGain(const SignalTy& window) {
	using Scalar = remove_complex_t<typename SignalTy::value_type>;
	return SumSquare(window) / Scalar(window.size());
}


//------------------------------------------------------------------------------
// List of window functions.
//------------------------------------------------------------------------------

template <mutable_signal_or_view_r SignalOut>
void HammingWindow(SignalOut&& out) {
	using R = scalar_type_t<std::decay_t<SignalOut>>;
	using U = remove_complex_t<R>;

	LinSpace(out, U(0), U(2) * pi_v<U>, true);
	Cos(out, out);
	out *= U(-0.46);
	out += U(0.54);
}

template <mutable_signal_or_view_r SignalOut>
void FlatTopWindow(SignalOut&& out) {
	using R = scalar_type_t<std::decay_t<SignalOut>>;
	using U = remove_complex_t<R>;

	U c0 = U(0.21557895);
	U c1 = U(-0.41663158);
	U c2 = U(0.277263158);
	U c3 = U(-0.083578947);
	U c4 = U(0.006947368);

	const U N = U(out.size());
	U preSize1 = U(2) * pi_v<U> / (N - U(1));
	U preSize2 = U(4) * pi_v<U> / (N - U(1));
	U preSize3 = U(6) * pi_v<U> / (N - U(1));
	U preSize4 = U(8) * pi_v<U> / (N - U(1));

	LinSpace(out, U(0), U(out.size() - 1), true);
	std::for_each(out.begin(), out.end(), [&](R& k) {
		const U kreal = std::real(k);
		k = c0
			+ c1 * std::cos(preSize1 * kreal)
			+ c2 * std::cos(preSize2 * kreal)
			+ c3 * std::cos(preSize3 * kreal)
			+ c4 * std::cos(preSize4 * kreal);
	});
}

template <mutable_signal_or_view_r SignalOut>
void RectangularWindow(SignalOut&& out) {
	using R = scalar_type_t<std::decay_t<SignalOut>>;
	using U = remove_complex_t<R>;
	std::fill(out.begin(), out.end(), R(U(1.0)));
}

template <mutable_signal_or_view_r SignalOut>
void TriangularWindow(SignalOut&& out) {
	using R = scalar_type_t<std::decay_t<SignalOut>>;
	using U = remove_complex_t<R>;
	LinSpace(out, U(0), U(2), true);
	out -= U(1);
	Abs(out, out);
	out *= U(-1);
	out += U(1);
}

template <mutable_signal_or_view_r SignalOut>
void BlackmanWindow(SignalOut&& out) {
	using R = scalar_type_t<std::decay_t<SignalOut>>;
	using U = remove_complex_t<R>;
	LinSpace(out, U(0), U(2) * pi_v<U>, true);
	std::for_each(out.begin(), out.end(), [&](R& k) {
		const U kreal = std::real(k);
		k = U(0.42) - U(0.5) * std::cos(kreal) + U(0.08) * std::cos(2 * kreal);
	});
}

template <mutable_signal_or_view_r SignalOut>
void BlackmanHarrisWindow(SignalOut&& out) {
	using R = scalar_type_t<std::decay_t<SignalOut>>;
	using U = remove_complex_t<R>;
	LinSpace(out, U(0), U(2) * pi_v<U>, true);
	std::for_each(out.begin(), out.end(), [&](R& k) {
		const U kreal = std::real(k);
		k = U(0.35875) - U(0.48829) * std::cos(kreal) + U(0.14128) * std::cos(2 * kreal) + U(-0.01168) * std::cos(3 * kreal);
	});
}

template <mutable_signal_or_view_r SignalOut, class Number>
void GaussianWindow(SignalOut&& out, Number sigma = 1.f) {
	using Scalar = scalar_type_t<std::decay_t<SignalOut>>;
	using Real = remove_complex_t<Scalar>;
	const auto N = Real(out.size());
	const auto M = (N - Real(1)) / Real(2);
	LinSpace(out, -M, M, true);
	out *= Real(1) / (Real(sigma) * M);
	Multiply(out, out, out);
	out *= Real(-0.5);
	Exp(out, out);
}

template <mutable_signal_or_view_r SignalOut, class Number>
void KaiserWindow(SignalOut&& out, Number alpha) {
	using Scalar = scalar_type_t<std::decay_t<SignalOut>>;
	using Real = remove_complex_t<Scalar>;
	LinSpace(out, -Real(1), Real(1), true);
	std::for_each(out.begin(), out.end(), [&](Scalar& k) {
		const Real kreal = std::real(k);
		const Real piAlpha = pi_v<Real> * Real(alpha);
		const Real arg = std::sqrt(std::max(Real(0), Real(1) - kreal * kreal));
		k = Real(std::cyl_bessel_i(Real(0), piAlpha * arg)) / Real(std::cyl_bessel_i(Real(0), Real(piAlpha)));
	});
}

template <mutable_signal_or_view_r SignalOut>
void LanczosWindow(SignalOut&& out) {
	using Scalar = scalar_type_t<std::decay_t<SignalOut>>;
	using Real = remove_complex_t<Scalar>;
	LinSpace(out, -pi_v<Real>, pi_v<Real>, true);
	std::for_each(out.begin(), out.end(), [&](Scalar& k) {
		const Real kreal = std::real(k);
		k = kreal != Real(0) ? std::sin(kreal) / kreal : Real(1);
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
	}

} // namespace impl

template <mutable_signal_or_view_r SignalOut, class Number>
	requires(!is_complex_v<typename std::decay_t<SignalOut>::value_type>)
void DolphChebyshevWindow(SignalOut&& out, Number attenuation) {
	using Scalar = typename std::decay_t<SignalOut>::value_type;
	using Real = remove_complex_t<Scalar>;

	const size_t M = out.size() - 1;
	const Real beta = std::cosh(Real(1) / M * std::acosh(Real(1) / Real(attenuation)));
	Spectrum<std::complex<Real>> spectrum(out.size() / 2 + 1);
	LinSpace(spectrum, Real(0), pi_v<Real> * (Real(spectrum.size() - 1) / Real(out.size())), true);
	std::for_each(spectrum.begin(), spectrum.end(), [M, beta](auto& k) {
		const auto i = std::complex<Real>(0, 1);
		const auto phase = std::exp(i * k * Real(M % 2));
		const auto amplitude = impl::ChebyshevPoly(M, beta * std::cos(std::real(k)));
		k = phase * amplitude;
	});

	Ifft(out, spectrum);
	FftShift(out, out);
	const Real normalization = kernels::TransformReduce(
		out.begin(), out.end(), Real(0),
		[](const auto& acc, const auto& v) { return kernels::math_functions::max(acc, v); },
		[](const auto& v) { return kernels::math_functions::abs(v); });
	out *= Real(1) / normalization;
}

template <mutable_signal_or_view_r SignalOut, class Number>
void DolphChebyshevWindow(SignalOut&& out, Number attenuation) {
	using R = typename std::decay_t<SignalOut>::value_type;
	using T = remove_complex_t<R>;
	constexpr auto domain = domain_v<std::decay_t<SignalOut>>;

	BasicSignal<T, domain> outReal(out.size());
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


/// <summary> The windows functions as functors. </summary>
/// <remarks> These functors can be passed to some functions, like filter descriptors,
///		to make windowing simpler. </remarks>
namespace windows {
	struct Hamming {
		template <mutable_signal_or_view_r SignalOut>
		auto operator()(SignalOut&& out) const {
			return HammingWindow(out);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return HammingWindow<T, Domain>(length);
		}
	} inline constexpr hamming;

	struct Rectangular {
		template <mutable_signal_or_view_r SignalOut>
		auto operator()(SignalOut&& out) const {
			return RectangularWindow(out);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return RectangularWindow<T, Domain>(length);
		}
	} inline constexpr rectangular;

	struct Flattop {
		template <mutable_signal_or_view_r SignalOut>
		auto operator()(SignalOut&& out) const {
			return FlatTopWindow(out);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return FlatTopWindow<T, Domain>(length);
		}
	} inline constexpr flattop;

	struct Triangular {
		template <mutable_signal_or_view_r SignalOut>
		auto operator()(SignalOut&& out) const {
			return TriangularWindow(out);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return TriangularWindow<T, Domain>(length);
		}
	} inline constexpr triangular;

	struct Blackman {
		template <mutable_signal_or_view_r SignalOut>
		auto operator()(SignalOut&& out) const {
			return BlackmanWindow(out);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return BlackmanWindow<T, Domain>(length);
		}
	} inline constexpr blackman;

	struct BlackmanHarris {
		template <mutable_signal_or_view_r SignalOut>
		auto operator()(SignalOut&& out) const {
			return BlackmanHarrisWindow(out);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return BlackmanHarrisWindow<T, Domain>(length);
		}
	} inline constexpr blackmanHarris;

	struct Gaussian {
		template <mutable_signal_or_view_r SignalOut>
		auto operator()(SignalOut&& out) const {
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
	} inline constexpr gaussian;

	struct Kaiser {
		template <mutable_signal_or_view_r SignalOut>
		auto operator()(SignalOut&& out) const {
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
	} inline constexpr kaiser;

	struct Lanczos {
		template <mutable_signal_or_view_r SignalOut>
		auto operator()(SignalOut&& out) const {
			return LanczosWindow(out);
		}
		template <class T, eSignalDomain Domain = eSignalDomain::TIME>
		auto operator()(size_t length) const {
			return LanczosWindow<T, Domain>(length);
		}
	} inline constexpr lanczos;

	struct DolphChebyshev {
		template <mutable_signal_or_view_r SignalOut>
		auto operator()(SignalOut&& out) const {
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
	} inline constexpr dolphChebyshev;

} // namespace windows


} // namespace dspbb