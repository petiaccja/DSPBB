#pragma once

#include "FIR/BandTransforms.hpp"
#include "FIR/Descs.hpp"
#include "FIR/Filter.hpp"
#include "FIR/LeastSquares.hpp"
#include "FIR/Windowed.hpp"
#include "FilterUtility.hpp"

#include <algorithm>


namespace dspbb {

//------------------------------------------------------------------------------
// Window method
//------------------------------------------------------------------------------

// Lowpass
template <class SignalR, class ParamType, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR&& out, const impl::LowpassDesc<impl::FirMethodWindowed, ParamType, WindowType>& desc) {
	fir::KernelWindowedLowpass(out, desc.cutoff, desc.window);
}

// Highpass
template <class SignalR, class ParamType, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR&& out, const impl::HighpassDesc<impl::FirMethodWindowed, ParamType, WindowType>& desc) {
	FirFilter(out, Lowpass(WINDOWED).Cutoff(desc.cutoff).Window(desc.window));
	fir::ComplementaryResponse(out, out);
}

// Bandpass
template <class SignalR, class ParamType, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR&& out, const impl::BandpassDesc<impl::FirMethodWindowed, ParamType, WindowType>& desc) {
	const ParamType bandWidth = desc.upper - desc.lower;
	const ParamType bandCenter = (desc.upper + desc.lower) / ParamType(2);
	FirFilter(out, Lowpass(WINDOWED).Cutoff(bandWidth / ParamType(2)).Window(desc.window));
	fir::ShiftResponse(out, out, bandCenter);
}

// Bandstop
template <class SignalR, class ParamType, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR&& out, const impl::BandstopDesc<impl::FirMethodWindowed, ParamType, WindowType>& desc) {
	FirFilter(out, Bandpass(WINDOWED).Band(desc.lower, desc.upper).Window(desc.window));
	fir::ComplementaryResponse(out, out);
}

// Arbitrary
template <class SignalR, class ResponseFunc, class WindowType, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void FirFilter(SignalR&& out, const impl::ArbitraryDesc<impl::FirMethodWindowed, ResponseFunc, WindowType>& desc) {
	fir::KernelWindowedArbitrary(out, desc.responseFunc, desc.window);
}

//------------------------------------------------------------------------------
// Least-squares method
//------------------------------------------------------------------------------

namespace impl {

	template <class T>
	T Smoothstep(const T& x) {
		const T c = std::clamp(x, T(0), T(1));
		return T(3) * c * c - T(2) * c * c * c;
	}

	template <class T>
	T LerpParam(T x, T lower, T upper) {
		return (x - lower) / (upper - lower);
	}

	template <class F, class Desc>
	auto LeastSquaresSplitWeight(F frequency, const Desc& desc) {
		if (frequency <= static_cast<F>(desc.cutoffBegin)) {
			return static_cast<F>(desc.weightLow);
		}
		if (frequency <= static_cast<F>(desc.cutoffEnd)) {
			return static_cast<F>(desc.weightTransition);
		}
		return static_cast<F>(desc.weightHigh);
	}

	template <class F, class Desc>
	auto LeastSquaresBandWeight(F frequency, const Desc& desc) {
		if (frequency <= static_cast<F>(desc.lowerBegin)) {
			return static_cast<F>(desc.weightLow);
		}
		if (frequency <= static_cast<F>(desc.lowerEnd)) {
			return static_cast<F>(desc.weightTransition1);
		}
		if (frequency <= static_cast<F>(desc.upperBegin)) {
			return static_cast<F>(desc.weightMid);
		}
		if (frequency <= static_cast<F>(desc.upperEnd)) {
			return static_cast<F>(desc.weightTransition2);
		}
		return static_cast<F>(desc.weightHigh);
	}

	template <class ParamType>
	auto TranslateLeastSquares(const LowpassDesc<FirMethodLeastSquares, ParamType>& desc) {
		const auto response = [desc](auto f) {
			using F = std::decay_t<decltype(f)>;
			return Smoothstep(LerpParam(f, F(desc.cutoffEnd), F(desc.cutoffBegin)));
		};
		const auto weight = [desc](auto f) {
			return LeastSquaresSplitWeight(f, desc);
		};
		return std::make_tuple(response, weight);
	}

	template <class ParamType>
	auto TranslateLeastSquares(const HighpassDesc<FirMethodLeastSquares, ParamType>& desc) {
		const auto response = [desc](auto f) {
			using F = std::decay_t<decltype(f)>;
			return Smoothstep(LerpParam(f, F(desc.cutoffBegin), F(desc.cutoffEnd)));
		};
		const auto weight = [desc](auto f) {
			return LeastSquaresSplitWeight(f, desc);
		};
		return std::make_tuple(response, weight);
	}

	template <class ParamType>
	auto TranslateLeastSquares(const BandpassDesc<FirMethodLeastSquares, ParamType>& desc) {
		const ParamType fmid = (desc.lowerEnd + desc.upperBegin) / ParamType(2);
		const auto response = [desc, fmid](auto f) {
			using F = std::decay_t<decltype(f)>;
			return f < fmid ? Smoothstep(LerpParam(f, F(desc.lowerBegin), F(desc.lowerEnd)))
							: Smoothstep(LerpParam(f, F(desc.upperEnd), F(desc.upperBegin)));
		};
		const auto weight = [desc](auto f) {
			return LeastSquaresBandWeight(f, desc);
		};
		return std::make_tuple(response, weight);
	}

	template <class ParamType>
	auto TranslateLeastSquares(const BandstopDesc<FirMethodLeastSquares, ParamType>& desc) {
		const ParamType fmid = (desc.lowerEnd + desc.upperBegin) / ParamType(2);
		const auto response = [desc, fmid](auto f) {
			using F = std::decay_t<decltype(f)>;
			return f < fmid ? Smoothstep(LerpParam(f, F(desc.lowerEnd), F(desc.lowerBegin)))
							: Smoothstep(LerpParam(f, F(desc.upperBegin), F(desc.upperEnd)));
		};
		const auto weight = [desc](auto f) {
			return LeastSquaresBandWeight(f, desc);
		};
		return std::make_tuple(response, weight);
	}

	template <class ResponseFunc, class WeightFunc>
	auto TranslateLeastSquares(const ArbitraryDesc<FirMethodLeastSquares, ResponseFunc, WeightFunc>& desc) {
		return std::make_tuple(desc.responseFunc, desc.weightFunc);
	}

} // namespace impl

template <class SignalR, template <typename, typename...> class Desc, class... Params, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
auto FirFilter(SignalR&& out, const Desc<impl::FirMethodLeastSquares, Params...>& desc)
	-> decltype(void(impl::TranslateLeastSquares(desc))) {
	const auto [response, weight] = impl::TranslateLeastSquares(desc);
	fir::KernelLeastSquares(out, response, weight, desc.grid);
}


//------------------------------------------------------------------------------
// Hilbert
//------------------------------------------------------------------------------

namespace impl {

	template <class WindowType>
	auto TranslateHilbert2HalfbandDesc(const HilbertDesc<FirMethodWindowed, WindowType>& desc) {
		return Lowpass(WINDOWED).Cutoff(0.5f).Window(desc.window);
	}

	template <class ParamType>
	auto TranslateHilbert2HalfbandDesc(const HilbertDesc<FirMethodLeastSquares, ParamType>& desc) {
		const ParamType transitionBand = desc.transitionWidth;
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
	BasicSignal<T, Domain> out(taps);
	FirFilter(out, responseDesc);
	return out;
}


} // namespace dspbb