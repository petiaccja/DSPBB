#pragma once

#include "../../Primitives/Signal.hpp"
#include "../Windowing.hpp"

#include <type_traits>


namespace dspbb {

namespace impl {
	struct FirMethod {};
	struct FirMethodWindowed : FirMethod {};
	struct FirMethodLeastSquares : FirMethod {};


	template <class Method, class... Params>
	struct LowpassDesc;

	template <class Method, class... Params>
	struct HighpassDesc;

	template <class Method, class... Params>
	struct BandpassDesc;

	template <class Method, class... Params>
	struct BandstopDesc;

	template <class Method, class... Params>
	struct ArbitraryDesc;

	template <class Method, class... Params>
	struct HilbertDesc;

	inline auto DefaultWeight = [](auto f) { return decltype(f)(1); };
	inline auto DefaultResponse = [](auto f) { return decltype(f)(1); };


	//------------------------------------------------------------------------------
	// Windowed descs
	//------------------------------------------------------------------------------

	template <template <typename, typename...> class Desc, class ParamType, class WindowType>
	struct SplitDescWindowed {
		ParamType cutoff = ParamType(0.5);
		WindowType window;

		template <class NewParamType>
		[[nodiscard]] auto Cutoff(NewParamType cutoffNew) {
			return Desc<FirMethodWindowed, NewParamType, WindowType>{ { std::move(cutoffNew), std::move(window) } };
		}
		template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, Signal<float, TIME_DOMAIN>&>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) {
			return Desc<FirMethodWindowed, ParamType, NewWindowType>{ { std::move(cutoff), std::move(windowNew) } };
		}
		template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) {
			return Desc<FirMethodWindowed, ParamType, NewWindowType>{ { std::move(cutoff), std::move(windowNew) } };
		}
	};

	template <template <typename, typename...> class Desc, class ParamType, class WindowType>
	struct BandDescWindowed {
		ParamType low = ParamType(0.25);
		ParamType high = ParamType(0.75);
		WindowType window;

		template <class NewParamType>
		[[nodiscard]] auto Band(NewParamType lowNew, NewParamType highNew) {
			return Desc<FirMethodWindowed, NewParamType, WindowType>{ { std::move(lowNew), std::move(highNew), std::move(window) } };
		}
		template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, Signal<float, TIME_DOMAIN>&>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) {
			return Desc<FirMethodWindowed, ParamType, NewWindowType>{ { std::move(low), std::move(high), std::move(windowNew) } };
		}
		template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) {
			return Desc<FirMethodWindowed, ParamType, NewWindowType>{ { std::move(low), std::move(high), std::move(windowNew) } };
		}
	};

	template <class ResponseFunc, class WindowType>
	struct ArbitraryDesc<FirMethodWindowed, ResponseFunc, WindowType> {
		ResponseFunc responseFunc;
		WindowType window;

		template <class NewResponseFunc, std::enable_if_t<std::is_invocable_v<NewResponseFunc, float>, int> = 0>
		[[nodiscard]] auto Response(NewResponseFunc responseFuncNew) {
			return ArbitraryDesc<FirMethodWindowed, NewResponseFunc, WindowType>{ std::move(responseFuncNew), std::move(window) };
		}
		template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, Signal<float, TIME_DOMAIN>&>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) {
			return ArbitraryDesc<FirMethodWindowed, ResponseFunc, NewWindowType>{ std::move(responseFunc), std::move(windowNew) };
		}
		template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) {
			return ArbitraryDesc<FirMethodWindowed, ResponseFunc, NewWindowType>{ std::move(responseFunc), std::move(windowNew) };
		}
	};

	template <class WindowType>
	struct HilbertDesc<FirMethodWindowed, WindowType> {
		WindowType window;

		template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, Signal<float, TIME_DOMAIN>&>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) {
			return HilbertDesc<FirMethodWindowed, NewWindowType>{ std::move(windowNew) };
		}
		template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) {
			return HilbertDesc<FirMethodWindowed, NewWindowType>{ std::move(windowNew) };
		}
	};

	template <class T, class WindowType>
	struct LowpassDesc<FirMethodWindowed, T, WindowType>
		: SplitDescWindowed<LowpassDesc, T, WindowType> {};

	template <>
	struct LowpassDesc<FirMethodWindowed>
		: SplitDescWindowed<LowpassDesc, float, windows::Hamming> {};

	template <class T, class WindowType>
	struct HighpassDesc<FirMethodWindowed, T, WindowType>
		: SplitDescWindowed<HighpassDesc, T, WindowType> {};

	template <>
	struct HighpassDesc<FirMethodWindowed>
		: SplitDescWindowed<HighpassDesc, float, windows::Hamming> {};

	template <class T, class WindowType>
	struct BandpassDesc<FirMethodWindowed, T, WindowType>
		: BandDescWindowed<BandpassDesc, T, WindowType> {};

	template <>
	struct BandpassDesc<FirMethodWindowed>
		: BandDescWindowed<BandpassDesc, float, windows::Hamming> {};

	template <class T, class WindowType>
	struct BandstopDesc<FirMethodWindowed, T, WindowType>
		: BandDescWindowed<BandstopDesc, T, WindowType> {};

	template <>
	struct BandstopDesc<FirMethodWindowed>
		: BandDescWindowed<BandstopDesc, float, windows::Hamming> {};

	template <>
	struct ArbitraryDesc<FirMethodWindowed>
		: ArbitraryDesc<FirMethodWindowed, decltype(DefaultResponse), windows::Hamming> {
		ArbitraryDesc() : ArbitraryDesc<FirMethodWindowed, decltype(DefaultResponse), windows::Hamming>{ DefaultResponse, windows::hamming } {}
	};

	template <>
	struct HilbertDesc<FirMethodWindowed>
		: HilbertDesc<FirMethodWindowed, windows::Hamming> {
		HilbertDesc() : HilbertDesc<FirMethodWindowed, windows::Hamming>{ windows::hamming } {}
	};

	//------------------------------------------------------------------------------
	// Least squares descs
	//------------------------------------------------------------------------------

	template <template <typename, typename...> class Desc, class ParamType>
	struct SplitDescLeastSquares {
		ParamType cutoffBegin = ParamType(0.45);
		ParamType cutoffEnd = ParamType(0.55);
		ParamType weightLow = ParamType(1.0);
		ParamType weightTransition = ParamType(0.0);
		ParamType weightHigh = ParamType(1.0);

		[[nodiscard]] auto Cutoff(ParamType begin, ParamType end) {
			return Desc<FirMethodLeastSquares, ParamType>{ { begin, end, weightLow, weightTransition, weightHigh } };
		}
		[[nodiscard]] auto Weight(ParamType low, ParamType transition, ParamType high) {
			return Desc<FirMethodLeastSquares, ParamType>{ { cutoffBegin, cutoffEnd, low, transition, high } };
		}
	};

	template <template <typename, typename...> class Desc, class ParamType>
	struct BandDescLeastSquares {
		ParamType cutoffBegin1 = ParamType(0.2);
		ParamType cutoffEnd1 = ParamType(0.3);
		ParamType cutoffBegin2 = ParamType(0.7);
		ParamType cutoffEnd2 = ParamType(0.8);
		ParamType weightLow = ParamType(1.0);
		ParamType weightTransition1 = ParamType(0.0);
		ParamType weightMid = ParamType(1.0);
		ParamType weightTransition2 = ParamType(0.0);
		ParamType weightHigh = ParamType(1.0);

		[[nodiscard]] auto Band(ParamType begin1, ParamType end1, ParamType begin2, ParamType end2) {
			return Desc<FirMethodLeastSquares, ParamType>{ { begin1, end1, begin2, end2, weightLow, weightTransition1, weightMid, weightTransition2, weightHigh } };
		}
		[[nodiscard]] auto Weight(ParamType low, ParamType transition1, ParamType mid, ParamType transition2, ParamType high) {
			return Desc<FirMethodLeastSquares, ParamType>{ { cutoffBegin1, cutoffEnd1, cutoffBegin2, cutoffEnd2, low, transition1, mid, transition2, high } };
		}
	};

	template <class ParamType>
	struct HilbertDesc<FirMethodLeastSquares, ParamType> {
		ParamType transition = ParamType(1.0);

		[[nodiscard]] auto TransitionWidth(ParamType bandwidthNew) {
			return HilbertDesc<FirMethodLeastSquares, ParamType>{ bandwidthNew };
		}
	};

	template <template <typename, typename...> class Desc>
	struct SplitDescLeastSquares<Desc, void> {
		template <class ParamType>
		[[nodiscard]] auto Cutoff(ParamType begin, ParamType end) {
			return Desc<FirMethodLeastSquares, ParamType>{}.Cutoff(begin, end);
		}
		template <class ParamType>
		[[nodiscard]] auto Weight(ParamType low, ParamType transition, ParamType high) {
			return Desc<FirMethodLeastSquares, ParamType>{}.Weight(low, transition, high);
		}
	};

	template <template <typename, typename...> class Desc>
	struct BandDescLeastSquares<Desc, void> {
		template <class ParamType>
		[[nodiscard]] auto Band(ParamType begin1, ParamType end1, ParamType begin2, ParamType end2) {
			return Desc<FirMethodLeastSquares, ParamType>{}.Band(begin1, end1, begin2, end2);
		}
		template <class ParamType>
		[[nodiscard]] auto Weight(ParamType low, ParamType transition1, ParamType mid, ParamType transition2, ParamType high) {
			return Desc<FirMethodLeastSquares, ParamType>{}.Weight(low, transition1, mid, transition2, high);
		}
	};

	template <class ResponseFunc, class WeightFunc>
	struct ArbitraryDesc<FirMethodLeastSquares, ResponseFunc, WeightFunc> {
		ResponseFunc responseFunc;
		WeightFunc weightFunc;

		template <class NewResponseFunc, std::enable_if_t<std::is_invocable_v<NewResponseFunc, float>, int> = 0>
		[[nodiscard]] auto Response(NewResponseFunc responseFuncNew) {
			return ArbitraryDesc<FirMethodLeastSquares, NewResponseFunc, WeightFunc>{ std::move(responseFuncNew), std::move(weightFunc) };
		}
		template <class NewWeightFunc, std::enable_if_t<std::is_invocable_v<NewWeightFunc, float>, int> = 0>
		[[nodiscard]] auto Weight(NewWeightFunc weightFuncNew) {
			return ArbitraryDesc<FirMethodLeastSquares, ResponseFunc, NewWeightFunc>{ std::move(responseFunc), std::move(weightFuncNew) };
		}
	};


	template <class T>
	struct LowpassDesc<FirMethodLeastSquares, T>
		: SplitDescLeastSquares<LowpassDesc, T> {};

	template <>
	struct LowpassDesc<FirMethodLeastSquares>
		: SplitDescLeastSquares<LowpassDesc, void> {};

	template <class T>
	struct HighpassDesc<FirMethodLeastSquares, T>
		: SplitDescLeastSquares<HighpassDesc, T> {};

	template <>
	struct HighpassDesc<FirMethodLeastSquares>
		: SplitDescLeastSquares<HighpassDesc, void> {};

	template <class T>
	struct BandpassDesc<FirMethodLeastSquares, T>
		: BandDescLeastSquares<BandpassDesc, T> {};

	template <>
	struct BandpassDesc<FirMethodLeastSquares>
		: BandDescLeastSquares<BandpassDesc, void> {};

	template <class T>
	struct BandstopDesc<FirMethodLeastSquares, T>
		: BandDescLeastSquares<BandstopDesc, T> {};

	template <>
	struct BandstopDesc<FirMethodLeastSquares>
		: BandDescLeastSquares<BandstopDesc, void> {};

	template <>
	struct ArbitraryDesc<FirMethodLeastSquares>
		: ArbitraryDesc<FirMethodLeastSquares, decltype(DefaultResponse), decltype(DefaultWeight)> {
		ArbitraryDesc() : ArbitraryDesc<FirMethodLeastSquares, decltype(DefaultResponse), decltype(DefaultWeight)>{ DefaultResponse, DefaultWeight } {}
	};

	template <>
	struct HilbertDesc<FirMethodLeastSquares>
		: HilbertDesc<FirMethodLeastSquares, float> {};

} // namespace impl

//------------------------------------------------------------------------------
// Factory functions
//------------------------------------------------------------------------------

constexpr impl::FirMethodWindowed WINDOWED;
constexpr impl::FirMethodLeastSquares LEAST_SQUARES;

template <class Method, std::enable_if_t<std::is_base_of_v<impl::FirMethod, Method>, int> = 0>
auto Lowpass(Method) {
	return impl::LowpassDesc<Method>{};
}

template <class Method, std::enable_if_t<std::is_base_of_v<impl::FirMethod, Method>, int> = 0>
auto Highpass(Method) {
	return impl::HighpassDesc<Method>{};
}

template <class Method, std::enable_if_t<std::is_base_of_v<impl::FirMethod, Method>, int> = 0>
auto Bandpass(Method) {
	return impl::BandpassDesc<Method>{};
}

template <class Method, std::enable_if_t<std::is_base_of_v<impl::FirMethod, Method>, int> = 0>
auto Bandstop(Method) {
	return impl::BandstopDesc<Method>{};
}

template <class Method, std::enable_if_t<std::is_base_of_v<impl::FirMethod, Method>, int> = 0>
auto Arbitrary(Method) {
	return impl::ArbitraryDesc<Method>{};
}

template <class Method, std::enable_if_t<std::is_base_of_v<impl::FirMethod, Method>, int> = 0>
auto Hilbert(Method) {
	return impl::HilbertDesc<Method>{};
}


} // namespace dspbb
