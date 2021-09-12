#pragma once

#include "FIRDescs.hpp"
#include "FIRLeastSquares.hpp"
#include "FIRTransforms.hpp"
#include "FIRWindowed.hpp"
#include "FIRCommon.hpp"


namespace dspbb {

//------------------------------------------------------------------------------
// Window method
//------------------------------------------------------------------------------

// Lowpass
template <class SignalR, class ParamType, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const LowpassDesc<MethodTagWindowed, ParamType, WindowType>& responseDesc) {
	FirLowpassWin(out, responseDesc.cutoff, responseDesc.window);
}


// Highpass
template <class SignalR, class ParamType, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const HighpassDesc<MethodTagWindowed, ParamType, WindowType>& responseDesc) {
	FirFilter(out, Lowpass(WINDOWED).Cutoff(responseDesc.cutoff).Window(responseDesc.window));
	ComplementaryResponse(out, out);
}

// Bandpass
template <class SignalR, class ParamType, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const BandpassDesc<MethodTagWindowed, ParamType, WindowType>& responseDesc) {
	const ParamType bandWidth = responseDesc.high - responseDesc.low;
	const ParamType bandCenter = (responseDesc.high + responseDesc.low) / ParamType(2);
	FirFilter(out, Lowpass(WINDOWED).Cutoff(bandWidth / ParamType(2)).Window(responseDesc.window));
	ShiftResponse(out, out, bandCenter);
}

// Bandstop
template <class SignalR, class ParamType, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const BandstopDesc<MethodTagWindowed, ParamType, WindowType>& responseDesc) {
	FirFilter(out, Bandpass(WINDOWED).Band(responseDesc.low, responseDesc.high).Window(responseDesc.window));
	ComplementaryResponse(out, out);
}

// Arbitrary
template <class SignalR, class ResponseFunc, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const ArbitraryDesc<MethodTagWindowed, ResponseFunc, WindowType>& responseDesc) {
	FirArbitraryWin(out, responseDesc.responseFunc, responseDesc.window);
}

// Hilbert
template <class SignalR, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const HilbertDesc<MethodTagWindowed, WindowType>& desc) {
	if (out.Size() % 2 == 0) {
		const size_t halfbandSize = out.Size() * 2 - 1;
		SignalR halfband(halfbandSize);
		FirFilter(halfband, Lowpass(WINDOWED).Cutoff(0.5f).Window(desc.window));
		HalfbandToHilbertEven(out, halfband);
	}
	else {
		FirFilter(out, Lowpass(WINDOWED).Cutoff(0.5f).Window(desc.window));
		HalfbandToHilbertOdd(out, out);
	}
}

//------------------------------------------------------------------------------
// Least-squares method
//------------------------------------------------------------------------------




//------------------------------------------------------------------------------
// Out-of-place wrapper.
//------------------------------------------------------------------------------

template <class T, eSignalDomain Domain, class ResponseDesc>
auto FirFilter(size_t taps, ResponseDesc responseDesc) {
	Signal<T, Domain> out(taps);
	FirFilter(out, responseDesc);
	return out;
}


} // namespace dspbb