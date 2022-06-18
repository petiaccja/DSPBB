#pragma once

#include "../../Primitives/Signal.hpp"
#include "../FilterUtility.hpp"
#include "../Windowing.hpp"

#include <type_traits>


namespace dspbb {

namespace impl {

	struct {
		template <class T>
		auto operator()(T f) const { return T(1); }
	} inline const DefaultResponse{};
	struct {
		template <class T>
		auto operator()(T f) const { return T(1); }
	} inline const DefaultWeight{};


	namespace windowed {

		template <template <typename, typename...> class Desc, class ParamType, class WindowType>
		struct SplitDescWindowed {
			ParamType cutoff = ParamType(0.5);
			WindowType window;

			template <class NewParamType>
			[[nodiscard]] auto Cutoff(NewParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<NewParamType, WindowType>{ { std::move(cutoffNew), window } };
			}
			template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, BasicSignal<float, TIME_DOMAIN>&>, int> = 0>
			[[nodiscard]] auto Window(NewWindowType windowNew) const {
				return Desc<ParamType, NewWindowType>{ { cutoff, std::move(windowNew) } };
			}
			template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
			[[nodiscard]] auto Window(NewWindowType windowNew) const {
				return Desc<ParamType, NewWindowType>{ { cutoff, std::move(windowNew) } };
			}
		};

		template <class T, class WindowType>
		struct LowpassDesc : SplitDescWindowed<LowpassDesc, T, WindowType> {};


		template <class T, class WindowType>
		struct HighpassDesc : SplitDescWindowed<HighpassDesc, T, WindowType> {};

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
				return Desc<NewParamType, WindowType>{ { std::move(lowerNew), std::move(upperNew), window } };
			}
			template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, BasicSignal<float, TIME_DOMAIN>&>, int> = 0>
			[[nodiscard]] auto Window(NewWindowType windowNew) const {
				return Desc<ParamType, NewWindowType>{ { lower, upper, std::move(windowNew) } };
			}
			template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
			[[nodiscard]] auto Window(NewWindowType windowNew) const {
				return Desc<ParamType, NewWindowType>{ { lower, upper, std::move(windowNew) } };
			}
		};

		template <class T, class WindowType>
		struct BandpassDesc : BandDescWindowed<BandpassDesc, T, WindowType> {};


		template <class T, class WindowType>
		struct BandstopDesc : BandDescWindowed<BandstopDesc, T, WindowType> {};

		template <class ResponseFunc, class WindowType>
		struct ArbitraryDesc {
			ResponseFunc responseFunc{};
			WindowType window;

			template <class NewResponseFunc, std::enable_if_t<std::is_invocable_v<NewResponseFunc, float>, int> = 0>
			[[nodiscard]] auto Response(NewResponseFunc responseFuncNew) const {
				return ArbitraryDesc<NewResponseFunc, WindowType>{ std::move(responseFuncNew), window };
			}
			template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, BasicSignal<float, TIME_DOMAIN>&>, int> = 0>
			[[nodiscard]] auto Window(NewWindowType windowNew) const {
				return ArbitraryDesc<ResponseFunc, NewWindowType>{ responseFunc, std::move(windowNew) };
			}
			template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
			[[nodiscard]] auto Window(NewWindowType windowNew) const {
				return ArbitraryDesc<ResponseFunc, NewWindowType>{ responseFunc, std::move(windowNew) };
			}
		};

		template <class WindowType>
		struct HilbertDesc {
			WindowType window;

			template <class NewWindowType, std::enable_if_t<!is_signal_like_v<NewWindowType> && std::is_invocable_v<WindowType, BasicSignal<float, TIME_DOMAIN>&>, int> = 0>
			[[nodiscard]] auto Window(NewWindowType windowNew) const {
				return HilbertDesc<NewWindowType>{ std::move(windowNew) };
			}
			template <class NewWindowType, std::enable_if_t<is_signal_like_v<NewWindowType>, int> = 0>
			[[nodiscard]] auto Window(NewWindowType windowNew) const {
				return HilbertDesc<NewWindowType>{ std::move(windowNew) };
			}
		};

	} // namespace windowed


	namespace least_squares {

		template <template <typename> class Desc, class ParamType>
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
				return Desc<ParamType>{ { cutoffBeginNew, cutoffEndNew, weightLow, weightTransition, weightHigh, grid } };
			}
			[[nodiscard]] auto Weight(ParamType newLow, ParamType newTransition, ParamType newHigh) const {
				return Desc<ParamType>{ { cutoffBegin, cutoffEnd, newLow, newTransition, newHigh, grid } };
			}
			[[nodiscard]] auto Grid(size_t gridNew) const {
				return Desc<ParamType>{ { cutoffBegin, cutoffEnd, weightLow, weightTransition, weightHigh, gridNew } };
			}
		};

		template <class T>
		struct LowpassDesc : SplitDescLeastSquares<LowpassDesc, T> {};


		template <class T>
		struct HighpassDesc : SplitDescLeastSquares<HighpassDesc, T> {};


		template <template <typename> class Desc, class ParamType>
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
				return Desc<ParamType>{ { lowerBeginNew, lowerEndNew, upperBeginNew, upperEndNew, weightLow, weightTransition1, weightMid, weightTransition2, weightHigh, grid } };
			}
			[[nodiscard]] auto Weight(ParamType lowNew, ParamType transition1New, ParamType midNew, ParamType transition2New, ParamType highNew) const {
				return Desc<ParamType>{ { lowerBegin, lowerEnd, upperBegin, upperEnd, lowNew, transition1New, midNew, transition2New, highNew, grid } };
			}
			[[nodiscard]] auto Grid(size_t gridNew) const {
				return Desc<ParamType>{ { lowerBegin, lowerEnd, upperBegin, upperEnd, weightLow, weightTransition1, weightMid, weightTransition2, weightHigh, gridNew } };
			}
		};


		template <class T>
		struct BandpassDesc : BandDescLeastSquares<BandpassDesc, T> {};


		template <class T>
		struct BandstopDesc : BandDescLeastSquares<BandstopDesc, T> {};


		template <class ParamType>
		struct HilbertDesc {
			ParamType transitionWidth = ParamType(1.0);
			size_t grid = 0;

			[[nodiscard]] auto TransitionWidth(ParamType newTransitionWidth) const {
				return HilbertDesc<ParamType>{ newTransitionWidth, grid };
			}
			[[nodiscard]] auto Grid(size_t gridNew) const {
				return HilbertDesc<ParamType>{ transitionWidth, gridNew };
			}
		};


		template <class ResponseFunc, class WeightFunc>
		struct ArbitraryDesc {
			ResponseFunc responseFunc{};
			WeightFunc weightFunc{};
			size_t grid = 0;

			template <class NewResponseFunc, std::enable_if_t<std::is_invocable_v<NewResponseFunc, float>, int> = 0>
			[[nodiscard]] auto Response(NewResponseFunc responseFuncNew) const {
				return ArbitraryDesc<NewResponseFunc, WeightFunc>{ std::move(responseFuncNew), weightFunc, grid };
			}
			template <class NewWeightFunc, std::enable_if_t<std::is_invocable_v<NewWeightFunc, float>, int> = 0>
			[[nodiscard]] auto Weight(NewWeightFunc weightFuncNew) const {
				return ArbitraryDesc<ResponseFunc, NewWeightFunc>{ responseFunc, std::move(weightFuncNew), grid };
			}
			[[nodiscard]] auto Grid(size_t gridNew) const {
				return ArbitraryDesc<ResponseFunc, WeightFunc>{ responseFunc, weightFunc, gridNew };
			}
		};


		template <template <typename> class Desc>
		struct SplitDescLeastSquares<Desc, void> {
			template <class ParamType>
			[[nodiscard]] auto Cutoff(ParamType cutoffBeginNew, ParamType cutoffEndNew) const {
				return Desc<ParamType>{}.Cutoff(cutoffBeginNew, cutoffEndNew);
			}
			template <class ParamType>
			[[nodiscard]] auto Weight(ParamType newLow, ParamType newTransition, ParamType newHigh) const {
				return Desc<ParamType>{}.Weight(newLow, newTransition, newHigh);
			}
		};

		template <>
		struct LowpassDesc<void> : SplitDescLeastSquares<LowpassDesc, void> {};


		template <>
		struct HighpassDesc<void> : SplitDescLeastSquares<HighpassDesc, void> {};


		template <template <typename> class Desc>
		struct BandDescLeastSquares<Desc, void> {
			template <class ParamType>
			[[nodiscard]] auto Band(ParamType lowerBeginNew, ParamType lowerEndNew, ParamType upperBeginNew, ParamType upperEndNew) const {
				return Desc<ParamType>{}.Band(lowerBeginNew, lowerEndNew, upperBeginNew, upperEndNew);
			}
			template <class ParamType>
			[[nodiscard]] auto Weight(ParamType lowNew, ParamType transition1New, ParamType midNew, ParamType transition2New, ParamType highNew) const {
				return Desc<ParamType>{}.Weight(lowNew, transition1New, midNew, transition2New, highNew);
			}
		};


		template <>
		struct BandpassDesc<void> : BandDescLeastSquares<BandpassDesc, void> {};


		template <>
		struct BandstopDesc<void> : BandDescLeastSquares<BandstopDesc, void> {};


		template <>
		struct HilbertDesc<void> {
			template <class ParamType>
			[[nodiscard]] auto TransitionWidth(ParamType newTransitionWidth) const {
				return HilbertDesc<ParamType>{}.TransitionWidth(newTransitionWidth);
			}
		};


	} // namespace least_squares

} // namespace impl

//------------------------------------------------------------------------------
// Factory functions
//------------------------------------------------------------------------------


struct {
	struct {
		const impl::windowed::LowpassDesc<float, windows::Hamming> Windowed{};
		const impl::least_squares::LowpassDesc<void> LeastSquares{};
	} const Lowpass{};
	struct {
		const impl::windowed::HighpassDesc<float, windows::Hamming> Windowed{};
		const impl::least_squares::HighpassDesc<void> LeastSquares{};
	} const Highpass{};
	struct {
		const impl::windowed::BandpassDesc<float, windows::Hamming> Windowed{};
		const impl::least_squares::BandpassDesc<void> LeastSquares{};
	} const Bandpass{};
	struct {
		const impl::windowed::BandstopDesc<float, windows::Hamming> Windowed{};
		const impl::least_squares::BandstopDesc<void> LeastSquares{};
	} const Bandstop{};
	struct {
		const impl::windowed::HilbertDesc<windows::Hamming> Windowed{};
		const impl::least_squares::HilbertDesc<void> LeastSquares{};
	} const Hilbert{};
	struct {
		const impl::windowed::ArbitraryDesc<decltype(impl::DefaultResponse), windows::Hamming> Windowed{};
		const impl::least_squares::ArbitraryDesc<decltype(impl::DefaultResponse), decltype(impl::DefaultWeight)> LeastSquares{};
	} const Arbitrary{};
} inline const Fir;


} // namespace dspbb
