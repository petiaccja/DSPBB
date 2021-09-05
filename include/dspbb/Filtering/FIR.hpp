#pragma once

#include "FIRCommon.hpp"
#include "FIRLeastSquares.hpp"
#include "FIRTransforms.hpp"
#include "FIRWindowed.hpp"


namespace dspbb {

//------------------------------------------------------------------------------
// Window method
//------------------------------------------------------------------------------

template <class SignalR, class U, class Func, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const LowpassDesc<U>& responseDesc, const WindowMethodFuncDesc<Func>& methodDesc) {
	FirLowpassWin(out, responseDesc.cutoff, methodDesc.windowFunc);
}

template <class SignalR, class U, class W, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const LowpassDesc<U>& responseDesc, const WindowMethodCoeffDesc<W, signal_traits<SignalR>::domain>& methodDesc) {
	FirLowpassWin(out, responseDesc.cutoff, methodDesc.windowCoefficients);
}

template <class SignalR, class U, class Func, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const ArbitraryFuncDesc<U>& responseDesc, const WindowMethodFuncDesc<Func>& methodDesc) {
	FirArbitraryWin(out, responseDesc.response, methodDesc.windowFunc);
}

template <class SignalR, class U, class W, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const ArbitraryFuncDesc<U>& responseDesc, const WindowMethodCoeffDesc<W, signal_traits<SignalR>::domain>& methodDesc) {
	FirArbitraryWin(out, responseDesc.response, methodDesc.windowCoefficients);
}


//------------------------------------------------------------------------------
// Least-squares method
//------------------------------------------------------------------------------

template <class SignalR, class U, class Func, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const LowpassDesc<U>& responseDesc, const LeastSquaresMethodFuncDesc<Func>& methodDesc) {
}

template <class SignalR, class U, class W, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const LowpassDesc<U>& responseDesc, const LeastSquaresMethodParamDesc<W>& methodDesc) {
	using R = typename std::decay_t<SignalR>::value_type;
	using T = remove_complex_t<R>;

	const auto ignoreLow = std::max(T(0), T(responseDesc.cutoff) - T(methodDesc.transitionBandwidth) / T(2));
	const auto ignoreHigh = std::min(T(1), T(responseDesc.cutoff) + T(methodDesc.transitionBandwidth) / T(2));
	const auto weightFunction = [&](float normalizedFrequency) {
		if (normalizedFrequency <= ignoreLow) {
			return T(methodDesc.passbandWeight);
		}
		if (normalizedFrequency <= ignoreHigh) {
			return T(0);
		}
		return T(methodDesc.stopbandWeight);
	};
	const auto responseFunction = [&](float normalizedFrequency) {
		if (normalizedFrequency < responseDesc.cutoff) {
			return T(1);
		}
		return T(0);
	};

	FirLeastSquares(out, responseFunction, weightFunction);
}

template <class SignalR, class U, class Func, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const ArbitraryFuncDesc<U>& responseDesc, const LeastSquaresMethodFuncDesc<Func>& methodDesc) {
	FirLeastSquares(out, responseDesc.response, methodDesc.weightFunction);
}



//------------------------------------------------------------------------------
// Band transforms
//------------------------------------------------------------------------------

// Highpass

template <class SignalR, class U, class MethodDesc, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const HighpassDesc<U>& responseDesc, const MethodDesc& methodDesc) {
	FirFilter(out, Lowpass(responseDesc.cutoff), methodDesc);
	ComplementaryResponse(out, out);
}

// Bandpass

template <class SignalR, class U, class MethodDesc, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const BandpassDesc<U>& responseDesc, const MethodDesc& methodDesc) {
	U bandWidth = responseDesc.upper - responseDesc.lower;
	U bandCenter = (responseDesc.upper + responseDesc.lower) / 2;
	FirFilter(out, Lowpass(bandWidth / U(2)), methodDesc);
	ShiftResponse(out, out, bandCenter);
}

// Bandstop

template <class SignalR, class U, class MethodDesc, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const BandstopDesc<U>& responseDesc, const MethodDesc& methodDesc) {
	FirFilter(out, Bandpass(responseDesc.lower, responseDesc.upper), methodDesc);
	ComplementaryResponse(out, out);
}

// Hilbert

template <class SignalR, class MethodDesc, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const HilbertDesc&, const MethodDesc& methodDesc) {
	if (out.Size() % 2 == 0) {
		const size_t halfbandSize = out.Size() * 2 - 1;
		SignalR halfband(halfbandSize);
		FirFilter(halfband, Lowpass(0.5f), methodDesc);
		HalfbandToHilbertEven(out, halfband);
	}
	else {
		FirFilter(out, Lowpass(0.5f), methodDesc);
		HalfbandToHilbertOdd(out, out);
	}
}

//------------------------------------------------------------------------------
// Out-of-place wrapper.
//------------------------------------------------------------------------------

template <class T, eSignalDomain Domain, class ResponseDesc, class MethodDesc>
auto FirFilter(size_t taps, ResponseDesc responseDesc, MethodDesc methodDesc) {
	Signal<T, Domain> out(taps);
	FirFilter(out, responseDesc, methodDesc);
	return out;
}


} // namespace dspbb