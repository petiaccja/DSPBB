#pragma once

#include "FIR/BandTransforms.hpp"
#include "FIR/Descs.hpp"
#include "FIR/LeastSquares.hpp"
#include "FIR/Utility.hpp"
#include "FIR/Windowed.hpp"

#include <algorithm>


namespace dspbb {

//------------------------------------------------------------------------------
// Window method
//------------------------------------------------------------------------------

// Lowpass
template <class SignalR, class ParamType, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR&& out, const impl::LowpassDesc<impl::MethodTagWindowed, ParamType, WindowType>& desc) {
	fir::KernelWindowedLowpass(out, desc.cutoff, desc.window);
}

// Highpass
template <class SignalR, class ParamType, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR&& out, const impl::HighpassDesc<impl::MethodTagWindowed, ParamType, WindowType>& desc) {
	FirFilter(out, Lowpass(WINDOWED).Cutoff(desc.cutoff).Window(desc.window));
	fir::ComplementaryResponse(out, out);
}

// Bandpass
template <class SignalR, class ParamType, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR&& out, const impl::BandpassDesc<impl::MethodTagWindowed, ParamType, WindowType>& desc) {
	const ParamType bandWidth = desc.high - desc.low;
	const ParamType bandCenter = (desc.high + desc.low) / ParamType(2);
	FirFilter(out, Lowpass(WINDOWED).Cutoff(bandWidth / ParamType(2)).Window(desc.window));
	fir::ShiftResponse(out, out, bandCenter);
}

// Bandstop
template <class SignalR, class ParamType, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR&& out, const impl::BandstopDesc<impl::MethodTagWindowed, ParamType, WindowType>& desc) {
	FirFilter(out, Bandpass(WINDOWED).Band(desc.low, desc.high).Window(desc.window));
	fir::ComplementaryResponse(out, out);
}

// Arbitrary
template <class SignalR, class ResponseFunc, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR&& out, const impl::ArbitraryDesc<impl::MethodTagWindowed, ResponseFunc, WindowType>& desc) {
	fir::KernelWindowedArbitrary(out, desc.responseFunc, desc.window);
}

//------------------------------------------------------------------------------
// Least-squares method
//------------------------------------------------------------------------------

namespace impl {

	template <class T>
	T Smoothstep(const T& x) {
		const float c = std::clamp(x, T(0), T(1));
		return T(3) * c * c - T(2) * c * c * c;
	}

	template <class F, class Desc>
	auto LeastSquaresSplitWeight(F frequency, const Desc& desc) {
		if (frequency <= desc.cutoffBegin) {
			return desc.weightLow;
		}
		if (frequency <= desc.cutoffEnd) {
			return desc.weightTransition;
		}
		return desc.weightHigh;
	}

	template <class F, class Desc>
	auto LeastSquaresBandWeight(F frequency, const Desc& desc) {
		if (frequency <= desc.cutoffBegin1) {
			return desc.weightLow;
		}
		if (frequency <= desc.cutoffEnd1) {
			return desc.weightTransition1;
		}
		if (frequency <= desc.cutoffBegin2) {
			return desc.weightMid;
		}
		if (frequency <= desc.cutoffEnd2) {
			return desc.weightTransition2;
		}
		return desc.weightHigh;
	}

	template <class ParamType>
	auto TranslateLeastSquares(const impl::LowpassDesc<impl::MethodTagLeastSquares, ParamType>& desc) {
		const auto response = [desc](ParamType f) {
			return Smoothstep((f - desc.cutoffEnd) / (desc.cutoffBegin - desc.cutoffEnd));
		};
		const auto weight = [desc](ParamType f) {
			return LeastSquaresSplitWeight(f, desc);
		};
		return std::make_tuple(response, weight);
	}

	template <class ParamType>
	auto TranslateLeastSquares(const impl::HighpassDesc<impl::MethodTagLeastSquares, ParamType>& desc) {
		const auto response = [desc](ParamType f) {
			return Smoothstep((f - desc.cutoffBegin) / (desc.cutoffEnd - desc.cutoffBegin));
		};
		const auto weight = [desc](ParamType f) {
			return LeastSquaresSplitWeight(f, desc);
		};
		return std::make_tuple(response, weight);
	}

	template <class ParamType>
	auto TranslateLeastSquares(const impl::BandpassDesc<impl::MethodTagLeastSquares, ParamType>& desc) {
		const ParamType fmid = (desc.cutoffEnd1 + desc.cutoffBegin2) / ParamType(2);
		const auto response = [desc, fmid](ParamType f) {
			return f < fmid ? Smoothstep((f - desc.cutoffBegin1) / (desc.cutoffEnd1 - desc.cutoffBegin1)) :
								Smoothstep((f - desc.cutoffEnd2) / (desc.cutoffBegin2 - desc.cutoffEnd2));
		};
		const auto weight = [desc](ParamType f) {
			return LeastSquaresBandWeight(f, desc);
		};
		return std::make_tuple(response, weight);
	}

	template <class ParamType>
	auto TranslateLeastSquares(const impl::BandstopDesc<impl::MethodTagLeastSquares, ParamType>& desc) {
		const ParamType fmid = (desc.cutoffEnd1 + desc.cutoffBegin2) / ParamType(2);
		const auto response = [desc, fmid](ParamType f) {
			return f < fmid ? Smoothstep((f - desc.cutoffEnd1) / (desc.cutoffBegin1 - desc.cutoffEnd1)) :
								Smoothstep((f - desc.cutoffBegin2) / (desc.cutoffEnd2 - desc.cutoffBegin2));
		};
		const auto weight = [desc](ParamType f) {
			return LeastSquaresBandWeight(f, desc);
		};
		return std::make_tuple(response, weight);
	}

	template <class ResponseFunc, class WeightFunc>
	auto TranslateLeastSquares(const impl::ArbitraryDesc<impl::MethodTagLeastSquares, ResponseFunc, WeightFunc>& desc) {
		return std::make_tuple(desc.responseFunc, desc.weightFunc);
	}

} // namespace impl

template <class SignalR, template <typename, typename...> class Desc, class... Params>
auto FirFilter(SignalR&& out, const Desc<impl::MethodTagLeastSquares, Params...>& desc)
	-> decltype(void(impl::TranslateLeastSquares(desc))) {
	const auto [response, weight] = impl::TranslateLeastSquares(desc);
	fir::KernelLeastSquares(out, response, weight);
}


//------------------------------------------------------------------------------
// Hilbert
//------------------------------------------------------------------------------

namespace impl {

	template <class WindowType>
	auto TranslateHilbert2HalfbandDesc(const impl::HilbertDesc<impl::MethodTagWindowed, WindowType>& desc) {
		return Lowpass(WINDOWED).Cutoff(0.5f).Window(desc.window);
	}

	template <class ParamType>
	auto TranslateHilbert2HalfbandDesc(const impl::HilbertDesc<impl::MethodTagLeastSquares, ParamType>& desc) {
		const ParamType transitionBand = desc.transition;
		return Lowpass(LEAST_SQUARES).Cutoff(ParamType(0.5) - transitionBand, ParamType(0.5) + transitionBand);
	}

} // namespace impl

template <class SignalR, class Method, class... Params, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR&& out, const impl::HilbertDesc<Method, Params...>& desc) {
	const auto halfbandDesc = impl::TranslateHilbert2HalfbandDesc(desc);

	if (out.Size() % 2 == 0) {
		const size_t halfbandSize = out.Size() * 2 - 1;
		std::decay_t<SignalR> halfband(halfbandSize);
		FirFilter(halfband, halfbandDesc);
		fir::HalfbandToHilbertEven(out, halfband);
	}
	else {
		FirFilter(out, halfbandDesc);
		fir::HalfbandToHilbertOdd(out, out);
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