#pragma once

#include <type_traits>
#include "../Primitives/Signal.hpp"


namespace dspbb {

template <class T, class U>
T NormalizedFrequency(T frequency, U sampleRate) {
	return T(2) * frequency / T(sampleRate);
}

//------------------------------------------------------------------------------
// Response description
//------------------------------------------------------------------------------
template <class T>
struct LowpassDesc {
	T cutoff;
};

template <class T>
struct HighpassDesc {
	T cutoff;
};

template <class T>
struct BandpassDesc {
	T lower;
	T upper;
};

template <class T>
struct BandstopDesc {
	T lower;
	T upper;
};

template <class Func>
struct ArbitraryFuncDesc {
	Func response;
};

struct HilbertDesc {};


template <class T>
auto Lowpass(T cutoff) {
	return LowpassDesc<T>{ cutoff };
}

template <class T>
auto Highpass(T cutoff) {
	return HighpassDesc<T>{ cutoff };
}

template <class T>
auto Bandpass(T lower, T upper) {
	return BandpassDesc<T>{ lower, upper };
}

template <class T>
auto Bandstop(T lower, T upper) {
	return BandstopDesc<T>{ lower, upper };
}

template <class Func>
auto Arbitrary(Func response) {
	return ArbitraryFuncDesc<Func>{ std::move(response) };
}

inline auto Hilbert() {
	return HilbertDesc{};
}

//------------------------------------------------------------------------------
// Method description
//------------------------------------------------------------------------------

// Window method description

template <class Func>
struct WindowMethodFuncDesc {
	Func windowFunc;
};

template <class T, eSignalDomain Domain>
struct WindowMethodCoeffDesc {
	SignalView<const T, Domain> windowCoefficients;
};

template <class Func, std::enable_if_t<std::is_invocable_v<Func, SignalView<float, TIME_DOMAIN>>, int> = 0>
auto Windowed(Func windowFunction) {
	return WindowMethodFuncDesc<Func>{ std::move(windowFunction) };
}

template <class SignalT, std::enable_if_t<is_signal_like_v<SignalT>, int> = 0>
auto Windowed(const SignalT& windowCoefficients) {
	using T = typename signal_traits<SignalT>::type;
	constexpr auto Domain = signal_traits<SignalT>::domain;
	return WindowMethodCoeffDesc<T, Domain>{ AsConstView(windowCoefficients) };
}

// Least squares method description.

template <class T>
struct LeastSquaresMethodParamDesc {
	T transitionBandwidth;
	T passbandWeight;
	T stopbandWeight;
};

template <class WeightFunc>
struct LeastSquaresMethodFuncDesc {
	WeightFunc weightFunction;
};

template <class T = float>
auto LeastSquares(T transitionBandwidth = 0.0f, T passbandWeight = T(1.0), T stopbandWeight = T(1.0)) {
	return LeastSquaresMethodParamDesc<T>{ transitionBandwidth, passbandWeight, stopbandWeight };
}

template <class WeightFunc, class = std::result_of_t<WeightFunc(float)>>
auto LeastSquares(WeightFunc weightFunc) {
	return LeastSquaresMethodFuncDesc<WeightFunc>{ std::move(weightFunc) };
}



} // namespace dspbb