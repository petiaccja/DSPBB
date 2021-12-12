#pragma once

#include "../../Primitives/Signal.hpp"
#include "../FilterUtility.hpp"
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
		[[nodiscard]] auto Cutoff(NewParamType cutoffNew) const {
			impl::ThrowIfNotNormalized(cutoffNew);
			return Desc<FirMethodWindowed, NewParamType, WindowType>{ { std::move(cutoffNew), window } };
		}
		template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, BasicSignal<float, TIME_DOMAIN>&>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) const {
			return Desc<FirMethodWindowed, ParamType, NewWindowType>{ { cutoff, std::move(windowNew) } };
		}
		template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) const {
			return Desc<FirMethodWindowed, ParamType, NewWindowType>{ { cutoff, std::move(windowNew) } };
		}
	};

	template <template <typename, typename...> class Desc, class ParamType, class WindowType>
	struct BandDescWindowed {
		ParamType lower = ParamType(0.25);
		ParamType upper = ParamType(0.75);
		WindowType window;

		template <class NewParamType>
		[[nodiscard]] auto Band(NewParamType lowerNew, NewParamType upperNew) const {
			impl::ThrowIfNotNormalized(lowerNew);
			impl::ThrowIfNotNormalized(upperNew);
			impl::ThrowIfNotSorted(lowerNew, upperNew);
			return Desc<FirMethodWindowed, NewParamType, WindowType>{ { std::move(lowerNew), std::move(upperNew), window } };
		}
		template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, BasicSignal<float, TIME_DOMAIN>&>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) const {
			return Desc<FirMethodWindowed, ParamType, NewWindowType>{ { lower, upper, std::move(windowNew) } };
		}
		template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) const {
			return Desc<FirMethodWindowed, ParamType, NewWindowType>{ { lower, upper, std::move(windowNew) } };
		}
	};

	template <class ResponseFunc, class WindowType>
	struct ArbitraryDesc<FirMethodWindowed, ResponseFunc, WindowType> {
		ResponseFunc responseFunc;
		WindowType window;

		template <class NewResponseFunc, std::enable_if_t<std::is_invocable_v<NewResponseFunc, float>, int> = 0>
		[[nodiscard]] auto Response(NewResponseFunc responseFuncNew) const {
			return ArbitraryDesc<FirMethodWindowed, NewResponseFunc, WindowType>{ std::move(responseFuncNew), window };
		}
		template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, BasicSignal<float, TIME_DOMAIN>&>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) const {
			return ArbitraryDesc<FirMethodWindowed, ResponseFunc, NewWindowType>{ responseFunc, std::move(windowNew) };
		}
		template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) const {
			return ArbitraryDesc<FirMethodWindowed, ResponseFunc, NewWindowType>{ responseFunc, std::move(windowNew) };
		}
	};

	template <class WindowType>
	struct HilbertDesc<FirMethodWindowed, WindowType> {
		WindowType window;

		template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, BasicSignal<float, TIME_DOMAIN>&>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) const {
			return HilbertDesc<FirMethodWindowed, NewWindowType>{ std::move(windowNew) };
		}
		template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
		[[nodiscard]] auto Window(NewWindowType windowNew) const {
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
		size_t grid = 0;

		[[nodiscard]] auto Cutoff(ParamType cutoffBeginNew, ParamType cutoffEndNew) const {
			impl::ThrowIfNotNormalized(cutoffBeginNew);
			impl::ThrowIfNotNormalized(cutoffEndNew);
			impl::ThrowIfNotSorted(cutoffBeginNew, cutoffEndNew);
			return Desc<FirMethodLeastSquares, ParamType>{ { cutoffBeginNew, cutoffEndNew, weightLow, weightTransition, weightHigh, grid } };
		}
		[[nodiscard]] auto Weight(ParamType newLow, ParamType newTransition, ParamType newHigh) const {
			return Desc<FirMethodLeastSquares, ParamType>{ { cutoffBegin, cutoffEnd, newLow, newTransition, newHigh, grid } };
		}
		[[nodiscard]] auto Grid(size_t gridNew) const {
			return Desc<FirMethodLeastSquares, ParamType>{ { cutoffBegin, cutoffEnd, weightLow, weightTransition, weightHigh, gridNew } };
		}
	};

	template <template <typename, typename...> class Desc, class ParamType>
	struct BandDescLeastSquares {
		ParamType lowerBegin = ParamType(0.2);
		ParamType lowerEnd = ParamType(0.3);
		ParamType upperBegin = ParamType(0.7);
		ParamType upperEnd = ParamType(0.8);
		ParamType weightLow = ParamType(1.0);
		ParamType weightTransition1 = ParamType(0.0);
		ParamType weightMid = ParamType(1.0);
		ParamType weightTransition2 = ParamType(0.0);
		ParamType weightHigh = ParamType(1.0);
		size_t grid = 0;

		[[nodiscard]] auto Band(ParamType lowerBeginNew, ParamType lowerEndNew, ParamType upperBeginNew, ParamType upperEndNew) const {
			impl::ThrowIfNotNormalized(lowerBeginNew);
			impl::ThrowIfNotNormalized(lowerEndNew);
			impl::ThrowIfNotNormalized(upperBeginNew);
			impl::ThrowIfNotNormalized(upperEndNew);
			impl::ThrowIfNotSorted(lowerBeginNew, lowerEndNew, upperBeginNew, upperEndNew);
			return Desc<FirMethodLeastSquares, ParamType>{ { lowerBeginNew, lowerEndNew, upperBeginNew, upperEndNew, weightLow, weightTransition1, weightMid, weightTransition2, weightHigh, grid } };
		}
		[[nodiscard]] auto Weight(ParamType lowNew, ParamType transition1New, ParamType midNew, ParamType transition2New, ParamType highNew) const {
			return Desc<FirMethodLeastSquares, ParamType>{ { lowerBegin, lowerEnd, upperBegin, upperEnd, lowNew, transition1New, midNew, transition2New, highNew, grid } };
		}
		[[nodiscard]] auto Grid(size_t gridNew) const {
			return Desc<FirMethodLeastSquares, ParamType>{ { lowerBegin, lowerEnd, upperBegin, upperEnd, weightLow, weightTransition1, weightMid, weightTransition2, weightHigh, gridNew } };
		}
	};

	template <class ParamType>
	struct HilbertDesc<FirMethodLeastSquares, ParamType> {
		ParamType transitionWidth = ParamType(1.0);
		size_t grid = 0;

		[[nodiscard]] auto TransitionWidth(ParamType newTransitionWidth) const {
			return HilbertDesc<FirMethodLeastSquares, ParamType>{ newTransitionWidth, grid };
		}
		[[nodiscard]] auto Grid(size_t gridNew) const {
			return HilbertDesc<FirMethodLeastSquares, ParamType>{ transitionWidth, gridNew };
		}
	};

	template <template <typename, typename...> class Desc>
	struct SplitDescLeastSquares<Desc, void> {
		template <class ParamType>
		[[nodiscard]] auto Cutoff(ParamType begin, ParamType end) const {
			return Desc<FirMethodLeastSquares, ParamType>{}.Cutoff(begin, end);
		}
		template <class ParamType>
		[[nodiscard]] auto Weight(ParamType low, ParamType transition, ParamType high) const {
			return Desc<FirMethodLeastSquares, ParamType>{}.Weight(low, transition, high);
		}
		// Initializer with Grid is missing because then it cannot determine ParamType.
	};

	template <template <typename, typename...> class Desc>
	struct BandDescLeastSquares<Desc, void> {
		template <class ParamType>
		[[nodiscard]] auto Band(ParamType lowerBegin, ParamType lowerEnd, ParamType upperBegin, ParamType upperEnd) const {
			return Desc<FirMethodLeastSquares, ParamType>{}.Band(lowerBegin, lowerEnd, upperBegin, upperEnd);
		}
		template <class ParamType>
		[[nodiscard]] auto Weight(ParamType low, ParamType transition1, ParamType mid, ParamType transition2, ParamType high) const {
			return Desc<FirMethodLeastSquares, ParamType>{}.Weight(low, transition1, mid, transition2, high);
		}
	};

	template <>
	struct HilbertDesc<FirMethodLeastSquares, void> {
		template <class ParamType>
		[[nodiscard]] auto TransitionWidth(ParamType newTransitionWidth) const {
			return HilbertDesc<FirMethodLeastSquares, ParamType>{}.TransitionWidth(newTransitionWidth);
		}
	};

	template <class ResponseFunc, class WeightFunc>
	struct ArbitraryDesc<FirMethodLeastSquares, ResponseFunc, WeightFunc> {
		ResponseFunc responseFunc;
		WeightFunc weightFunc;
		size_t grid = 0;

		template <class NewResponseFunc, std::enable_if_t<std::is_invocable_v<NewResponseFunc, float>, int> = 0>
		[[nodiscard]] auto Response(NewResponseFunc responseFuncNew) const {
			return ArbitraryDesc<FirMethodLeastSquares, NewResponseFunc, WeightFunc>{ std::move(responseFuncNew), weightFunc, grid };
		}
		template <class NewWeightFunc, std::enable_if_t<std::is_invocable_v<NewWeightFunc, float>, int> = 0>
		[[nodiscard]] auto Weight(NewWeightFunc weightFuncNew) const {
			return ArbitraryDesc<FirMethodLeastSquares, ResponseFunc, NewWeightFunc>{ responseFunc, std::move(weightFuncNew), grid };
		}
		[[nodiscard]] auto Grid(size_t gridNew) const {
			return ArbitraryDesc<FirMethodLeastSquares, ResponseFunc, WeightFunc>{ responseFunc, weightFunc, gridNew };
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
