#pragma once

#include "../Generators/Spaces.hpp"
#include "../Math/Statistics.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "../Utility/Algorithm.hpp"
#include "../Utility/Numbers.hpp"
#include "../Utility/TypeTraits.hpp"

#include <cmath>


namespace dspbb {

//------------------------------------------------------------------------------
// Assess properties of windows.
//------------------------------------------------------------------------------
template <class T, eSignalDomain Domain>
T CoherentGain(SignalView<T, Domain> window) {
	return Sum(window) / remove_complex_t<T>(window.Size());
}

template <class T, eSignalDomain Domain>
T CoherentGain(const Signal<T, Domain>& window) {
	return CoherentGain(AsConstView(window));
}

template <class T, eSignalDomain Domain>
T EnergyGain(SignalView<T, Domain> window) {
	return SumSquare(window) / T(window.Size());
}

template <class T, eSignalDomain Domain>
T EnergyGain(const Signal<T, Domain>& window) {
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
		const U kreal = std::real(k); // We can do that safely for std::complex, it's a POD type.
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
		const U kreal = std::real(k); // We can do that safely for std::complex, it's a POD type.
		k = U(0.42) - U(0.5) * std::cos(kreal) + U(0.08) * std::cos(2 * kreal);
	});
}

template <class SignalR, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void BlackmanHarrisWindow(SignalR&& out) {
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using U = remove_complex_t<R>;
	LinSpace(out, U(0), U(2) * pi_v<U>, true);
	std::for_each(out.begin(), out.end(), [&](R& k) {
		const U kreal = std::real(k); // We can do that safely for std::complex, it's a POD type.
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


template <class T, eSignalDomain Domain = eSignalDomain::TIME>
Signal<T, Domain> HammingWindow(size_t length) {
	Signal<T, Domain> window(length);
	HammingWindow(AsView(window));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
Signal<T, Domain> FlatTopWindow(size_t length) {
	Signal<T, Domain> window(length);
	FlatTopWindow(AsView(window));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
Signal<T, Domain> RectangularWindow(size_t length) {
	Signal<T, Domain> window(length, T(1.0));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
Signal<T, Domain> TriangularWindow(size_t length) {
	Signal<T, Domain> window(length);
	TriangularWindow(AsView(window));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
Signal<T, Domain> BlackmanWindow(size_t length) {
	Signal<T, Domain> window(length);
	BlackmanWindow(AsView(window));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
Signal<T, Domain> BlackmanHarrisWindow(size_t length) {
	Signal<T, Domain> window(length);
	BlackmanHarrisWindow(AsView(window));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
Signal<T, Domain> GaussianWindow(size_t length, T sigma = 1) {
	Signal<T, Domain> window(length);
	GaussianWindow(AsView(window), sigma);
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
} // namespace windows


} // namespace dspbb