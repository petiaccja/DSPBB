#pragma once

#include "WindowFunctions.hpp"
#include "dspbb/Primitives/Signal.hpp"
#include "dspbb/Primitives/SignalView.hpp"

#include <type_traits>


namespace dspbb {

struct MethodTagWindowed {};
struct MethodTagLeastSquares {};

constexpr MethodTagWindowed WINDOWED;
constexpr MethodTagLeastSquares LEAST_SQUARES;

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
		return Desc<MethodTagWindowed, NewParamType, WindowType>{ { std::move(cutoffNew), std::move(window) } };
	}
	template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, Signal<float, TIME_DOMAIN>&>, int> = 0>
	[[nodiscard]] auto Window(NewWindowType windowNew) {
		return Desc<MethodTagWindowed, ParamType, NewWindowType>{ { std::move(cutoff), std::move(windowNew) } };
	}
	template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
	[[nodiscard]] auto Window(NewWindowType windowNew) {
		return Desc<MethodTagWindowed, ParamType, NewWindowType>{ { std::move(cutoff), std::move(windowNew) } };
	}
};

template <template <typename, typename...> class Desc, class ParamType, class WindowType>
struct BandDescWindowed {
	ParamType low = ParamType(0.25);
	ParamType high = ParamType(0.75);
	WindowType window;

	template <class NewParamType>
	[[nodiscard]] auto Band(NewParamType lowNew, NewParamType highNew) {
		return Desc<MethodTagWindowed, NewParamType, WindowType>{ { std::move(lowNew), std::move(highNew), std::move(window) } };
	}
	template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, Signal<float, TIME_DOMAIN>&>, int> = 0>
	[[nodiscard]] auto Window(NewWindowType windowNew) {
		return Desc<MethodTagWindowed, ParamType, NewWindowType>{ { std::move(low), std::move(high), std::move(windowNew) } };
	}
	template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
	[[nodiscard]] auto Window(NewWindowType windowNew) {
		return Desc<MethodTagWindowed, ParamType, NewWindowType>{ { std::move(low), std::move(high), std::move(windowNew) } };
	}
};

template <class ResponseFunc, class WindowType>
struct ArbitraryDesc<MethodTagWindowed, ResponseFunc, WindowType> {
	ResponseFunc responseFunc;
	WindowType window;

	template <class NewResponseFunc, std::enable_if_t<std::is_invocable_v<NewResponseFunc, float>, int> = 0>
	[[nodiscard]] auto Response(NewResponseFunc responseFuncNew) {
		return ArbitraryDesc<MethodTagWindowed, NewResponseFunc, WindowType>{ std::move(responseFuncNew), std::move(window) };
	}
	template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, Signal<float, TIME_DOMAIN>&>, int> = 0>
	[[nodiscard]] auto Window(NewWindowType windowNew) {
		return ArbitraryDesc<MethodTagWindowed, ResponseFunc, NewWindowType>{ std::move(responseFunc), std::move(windowNew) };
	}
	template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
	[[nodiscard]] auto Window(NewWindowType windowNew) {
		return ArbitraryDesc<MethodTagWindowed, ResponseFunc, NewWindowType>{ std::move(responseFunc), std::move(windowNew) };
	}
};

template <class WindowType>
struct HilbertDesc<MethodTagWindowed, WindowType> {
	WindowType window;
	
	template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, Signal<float, TIME_DOMAIN>&>, int> = 0>
	[[nodiscard]] auto Window(NewWindowType windowNew) {
		return HilbertDesc<MethodTagWindowed, NewWindowType>{ std::move(windowNew) };
	}
	template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
	[[nodiscard]] auto Window(NewWindowType windowNew) {
		return HilbertDesc<MethodTagWindowed, NewWindowType>{ std::move(windowNew) };
	}
};

template <class T, class WindowType>
struct LowpassDesc<MethodTagWindowed, T, WindowType>
	: SplitDescWindowed<LowpassDesc, T, WindowType> {};

template <>
struct LowpassDesc<MethodTagWindowed>
	: SplitDescWindowed<LowpassDesc, float, windows::Hamming> {};

template <class T, class WindowType>
struct HighpassDesc<MethodTagWindowed, T, WindowType>
	: SplitDescWindowed<HighpassDesc, T, WindowType> {};

template <>
struct HighpassDesc<MethodTagWindowed>
	: SplitDescWindowed<HighpassDesc, float, windows::Hamming> {};

template <class T, class WindowType>
struct BandpassDesc<MethodTagWindowed, T, WindowType>
	: BandDescWindowed<BandpassDesc, T, WindowType> {};

template <>
struct BandpassDesc<MethodTagWindowed>
	: BandDescWindowed<BandpassDesc, float, windows::Hamming> {};

template <class T, class WindowType>
struct BandstopDesc<MethodTagWindowed, T, WindowType>
	: BandDescWindowed<BandstopDesc, T, WindowType> {};

template <>
struct BandstopDesc<MethodTagWindowed>
	: BandDescWindowed<BandstopDesc, float, windows::Hamming> {};

template <>
struct ArbitraryDesc<MethodTagWindowed>
	: ArbitraryDesc<MethodTagWindowed, decltype(DefaultResponse), windows::Hamming> {
	ArbitraryDesc() : ArbitraryDesc<MethodTagWindowed, decltype(DefaultResponse), windows::Hamming>{ DefaultResponse, windows::hamming } {}
};

template <>
struct HilbertDesc<MethodTagWindowed>
	: HilbertDesc<MethodTagWindowed, windows::Hamming> {
	HilbertDesc() : HilbertDesc<MethodTagWindowed, windows::Hamming>{ windows::hamming } {}
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
	bool smooth = false;

	[[nodiscard]] auto Cutoff(ParamType begin, ParamType end) {
		return Desc<MethodTagLeastSquares, ParamType>{ { begin, end, weightLow, weightTransition, weightHigh, smooth } };
	}
	[[nodiscard]] auto Weight(ParamType low, ParamType transition, ParamType high) {
		return Desc<MethodTagLeastSquares, ParamType>{ { cutoffBegin, cutoffEnd, low, transition, high, smooth } };
	}
	[[nodiscard]] auto Smooth(bool enable) {
		return Desc<MethodTagLeastSquares, ParamType>{ { cutoffBegin, cutoffEnd, weightLow, weightTransition, weightHigh, enable } };
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
	bool smooth1 = false;
	bool smooth2 = false;

	[[nodiscard]] auto Band(ParamType begin1, ParamType end1, ParamType begin2, ParamType end2) {
		return Desc<MethodTagLeastSquares, ParamType>{ { begin1, end1, begin2, end2, weightLow, weightTransition1, weightMid, weightTransition2, weightHigh, smooth1, smooth2 } };
	}
	[[nodiscard]] auto Weight(ParamType low, ParamType transition1, ParamType mid, ParamType transition2, ParamType high) {
		return Desc<MethodTagLeastSquares, ParamType>{ { cutoffBegin1, cutoffEnd1, cutoffBegin2, cutoffEnd2, low, transition1, mid, transition2, high, smooth1, smooth2 } };
	}
	[[nodiscard]] auto Smooth(bool enable1, bool enable2) {
		return Desc<MethodTagLeastSquares, ParamType>{ { cutoffBegin1, cutoffEnd1, cutoffBegin2, cutoffEnd2, weightLow, weightTransition1, weightMid, weightTransition2, weightHigh, enable1, enable2 } };
	}
};

template <class ParamType>
struct HilbertDesc<MethodTagLeastSquares, ParamType> {
	ParamType transition = ParamType(1.0);
	
	[[nodiscard]] auto TransitionWidth(ParamType bandwidthNew) {
		return HilbertDesc<MethodTagLeastSquares, ParamType>{ bandwidthNew };
	}
};

template <template <typename, typename...> class Desc>
struct SplitDescLeastSquares<Desc, void> {
	template <class ParamType>
	[[nodiscard]] auto Cutoff(ParamType begin, ParamType end) {
		return Desc<MethodTagLeastSquares, ParamType>{}.Cutoff(begin, end);
	}
	template <class ParamType>
	[[nodiscard]] auto Weight(ParamType low, ParamType transition, ParamType high) {
		return Desc<MethodTagLeastSquares, ParamType>{}.Weight(low, transition, high);
	}
	template <class ParamType>
	[[nodiscard]] auto Smooth(bool enable) {
		return Desc<MethodTagLeastSquares, ParamType>{}.Smooth(enable);
	}
};

template <template <typename, typename...> class Desc>
struct BandDescLeastSquares<Desc, void> {
	template <class ParamType>
	[[nodiscard]] auto Band(ParamType begin1, ParamType end1, ParamType begin2, ParamType end2) {
		return Desc<MethodTagLeastSquares, ParamType>{}.Band(begin1, end1, begin2, end2);
	}
	template <class ParamType>
	[[nodiscard]] auto Weight(ParamType low, ParamType transition1, ParamType mid, ParamType transition2, ParamType high) {
		return Desc<MethodTagLeastSquares, ParamType>{}.Weight(low, transition1, mid, transition2, high);
	}
	template <class ParamType>
	[[nodiscard]] auto Smooth(bool enable1, bool enable2) {
		return Desc<MethodTagLeastSquares, ParamType>{}.Smooth(enable1, enable2);
	}
};

template <class ResponseFunc, class WeightFunc>
struct ArbitraryDesc<MethodTagLeastSquares, ResponseFunc, WeightFunc> {
	ResponseFunc responseFunc;
	WeightFunc weightFunc;

	template <class NewResponseFunc, std::enable_if_t<std::is_invocable_v<NewResponseFunc, float>, int> = 0>
	[[nodiscard]] auto Response(NewResponseFunc responseFuncNew) {
		return ArbitraryDesc<MethodTagLeastSquares, NewResponseFunc, WeightFunc>{ std::move(responseFuncNew), std::move(weightFunc) };
	}
	template <class NewWeightFunc, std::enable_if_t<std::is_invocable_v<NewWeightFunc, float>, int> = 0>
	[[nodiscard]] auto Weight(NewWeightFunc weightFuncNew) {
		return ArbitraryDesc<MethodTagLeastSquares, ResponseFunc, NewWeightFunc>{ std::move(responseFunc), std::move(weightFuncNew) };
	}
};


template <class T>
struct LowpassDesc<MethodTagLeastSquares, T>
	: SplitDescLeastSquares<LowpassDesc, T> {};

template <>
struct LowpassDesc<MethodTagLeastSquares>
	: SplitDescLeastSquares<LowpassDesc, void> {};

template <class T>
struct HighpassDesc<MethodTagLeastSquares, T>
	: SplitDescLeastSquares<HighpassDesc, T> {};

template <>
struct HighpassDesc<MethodTagLeastSquares>
	: SplitDescLeastSquares<HighpassDesc, void> {};

template <class T>
struct BandpassDesc<MethodTagLeastSquares, T>
	: BandDescLeastSquares<BandpassDesc, T> {};

template <>
struct BandpassDesc<MethodTagLeastSquares>
	: BandDescLeastSquares<BandpassDesc, void> {};

template <class T>
struct BandstopDesc<MethodTagLeastSquares, T>
	: BandDescLeastSquares<BandstopDesc, T> {};

template <>
struct BandstopDesc<MethodTagLeastSquares>
	: BandDescLeastSquares<BandstopDesc, void> {};

template <>
struct ArbitraryDesc<MethodTagLeastSquares>
	: ArbitraryDesc<MethodTagLeastSquares, decltype(DefaultResponse), decltype(DefaultWeight)> {
	ArbitraryDesc() : ArbitraryDesc<MethodTagLeastSquares, decltype(DefaultResponse), decltype(DefaultWeight)>{ DefaultResponse, DefaultWeight } {}
};

template <>
struct HilbertDesc<MethodTagLeastSquares>
	: HilbertDesc<MethodTagLeastSquares, float> {};

//------------------------------------------------------------------------------
// Factory functions
//------------------------------------------------------------------------------

template <class Method>
auto Lowpass(Method) {
	return LowpassDesc<Method>{};
}

template <class Method>
auto Highpass(Method) {
	return HighpassDesc<Method>{};
}

template <class Method>
auto Bandpass(Method) {
	return BandpassDesc<Method>{};
}

template <class Method>
auto Bandstop(Method) {
	return BandstopDesc<Method>{};
}

template <class Method>
auto Arbitrary(Method) {
	return ArbitraryDesc<Method>{};
}

template <class Method>
auto Hilbert(Method) {
	return HilbertDesc<Method>{};
}


} // namespace dspbb
