#pragma once


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

template <class Func>
struct WindowMethodFuncDesc {
	Func windowFunc;
};

template <class T, eSignalDomain Domain>
struct WindowMethodCoeffDesc {
	SignalView<const T, Domain> windowCoefficients;
};

template <class Func, class = std::result_of_t<Func(SignalView<float, TIME_DOMAIN>)>>
auto Windowed(Func windowFunction) {
	return WindowMethodFuncDesc<Func>{ std::move(windowFunction) };
}

template <class SignalT, std::enable_if_t<is_signal_like_v<SignalT>, int> = 0>
auto Windowed(const SignalT& windowCoefficients) {
	using T = typename signal_traits<SignalT>::type;
	constexpr auto Domain = signal_traits<SignalT>::domain;
	return WindowMethodCoeffDesc<T, Domain>{ AsConstView(windowCoefficients) };
}

} // namespace dspbb