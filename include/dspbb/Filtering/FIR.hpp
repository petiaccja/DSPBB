#pragma once

#include "FIRCommon.hpp"
#include "FIRDescs.hpp"
#include "FIRLeastSquares.hpp"
#include "FIRTransforms.hpp"
#include "FIRWindowed.hpp"

#include <algorithm>


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

//------------------------------------------------------------------------------
// Least-squares method
//------------------------------------------------------------------------------

template <class T>
T Smoothstep(const T& x) {
	const float c = std::clamp(x, T(0), T(1));
	return T(3) * c * c - T(2) * c * c * c;
}

template <class ParamType>
auto TranslateLeastSquares(const LowpassDesc<MethodTagLeastSquares, ParamType> desc) {
	const auto response = [=](ParamType f) {
		return Smoothstep((f - desc.cutoffEnd) / (desc.cutoffBegin - desc.cutoffEnd));
	};
	const auto weight = [=](ParamType f) {
		return f <= desc.cutoffBegin ? desc.weightLow :
			   f <= desc.cutoffEnd	 ? desc.weightTransition :
										 desc.weightHigh;
	};
	return std::make_tuple(response, weight);
}

template <class ParamType>
auto TranslateLeastSquares(const HighpassDesc<MethodTagLeastSquares, ParamType> desc) {
	const auto response = [=](ParamType f) {
		return Smoothstep((f - desc.cutoffBegin) / (desc.cutoffEnd - desc.cutoffBegin));
	};
	const auto weight = [=](ParamType f) {
		return f <= desc.cutoffBegin ? desc.weightLow :
			   f <= desc.cutoffEnd	 ? desc.weightTransition :
										 desc.weightHigh;
	};
	return std::make_tuple(response, weight);
}

template <class ParamType>
auto TranslateLeastSquares(const BandpassDesc<MethodTagLeastSquares, ParamType> desc) {
	const ParamType fmid = (desc.cutoffEnd1 + desc.cutoffBegin2) / ParamType(2);
	const auto response = [=](ParamType f) {
		return f < fmid ? Smoothstep((f - desc.cutoffBegin1) / (desc.cutoffEnd1 - desc.cutoffBegin1)) :
							Smoothstep((f - desc.cutoffEnd2) / (desc.cutoffBegin2 - desc.cutoffEnd2));
	};
	const auto weight = [=](ParamType f) {
		return f <= desc.cutoffBegin1 ? desc.weightLow :
			   f <= desc.cutoffEnd1	  ? desc.weightTransition1 :
			   f <= desc.cutoffBegin2 ? desc.weightMid :
			   f <= desc.cutoffEnd2	  ? desc.weightTransition2 :
										  desc.weightHigh;
	};
	return std::make_tuple(response, weight);
}

template <class ParamType>
auto TranslateLeastSquares(const BandstopDesc<MethodTagLeastSquares, ParamType> desc) {
	const ParamType fmid = (desc.cutoffEnd1 + desc.cutoffBegin2) / ParamType(2);
	const auto response = [=](ParamType f) {
		return f < fmid ? Smoothstep((f - desc.cutoffEnd1) / (desc.cutoffBegin1 - desc.cutoffEnd1)) :
							Smoothstep((f - desc.cutoffBegin2) / (desc.cutoffEnd2 - desc.cutoffBegin2));
	};
	const auto weight = [=](ParamType f) {
		return f <= desc.cutoffBegin1 ? desc.weightLow :
			   f <= desc.cutoffEnd1	  ? desc.weightTransition1 :
			   f <= desc.cutoffBegin2 ? desc.weightMid :
			   f <= desc.cutoffEnd2	  ? desc.weightTransition2 :
										  desc.weightHigh;
	};
	return std::make_tuple(response, weight);
}

template <class ResponseFunc, class WeightFunc>
auto TranslateLeastSquares(const ArbitraryDesc<MethodTagLeastSquares, ResponseFunc, WeightFunc> desc) {
	return std::make_tuple(desc.responseFunc, desc.weightFunc);
}

template <class SignalR, template <typename, typename...> class Desc, class... Params>
auto FirFilter(SignalR&& out, const Desc<MethodTagLeastSquares, Params...>& desc)
	-> decltype(void(TranslateLeastSquares(desc))) {
	const auto [response, weight] = TranslateLeastSquares(desc);
	FirLeastSquares(out, response, weight);
}


//------------------------------------------------------------------------------
// Hilbert
//------------------------------------------------------------------------------

template <class WindowType>
auto TranslateHilbert2HalfbandDesc(const HilbertDesc<MethodTagWindowed, WindowType>& desc) {
	return Lowpass(WINDOWED).Cutoff(0.5f).Window(desc.window);
}

template <class ParamType>
auto TranslateHilbert2HalfbandDesc(const HilbertDesc<MethodTagLeastSquares, ParamType>& desc) {
	const ParamType transitionBand = (ParamType(1) - desc.bandwidth) / ParamType(2);
	return Lowpass(LEAST_SQUARES).Cutoff(ParamType(0.5) - transitionBand, ParamType(0.5) + transitionBand);
}

template <class SignalR, class Method, class... Params, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR& out, const HilbertDesc<Method, Params...>& desc) {
	const auto halfbandDesc = TranslateHilbert2HalfbandDesc(desc);

	if (out.Size() % 2 == 0) {
		const size_t halfbandSize = out.Size() * 2 - 1;
		SignalR halfband(halfbandSize);
		FirFilter(halfband, halfbandDesc);
		HalfbandToHilbertEven(out, halfband);
	}
	else {
		FirFilter(out, halfbandDesc);
		HalfbandToHilbertOdd(out, out);
	}
}

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