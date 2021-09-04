#pragma once

#include "FIRWindowed.hpp"
#include "FIRTransforms.hpp"
#include "FIRLeastSquares.hpp"
#include "FIRCommon.hpp"


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


}