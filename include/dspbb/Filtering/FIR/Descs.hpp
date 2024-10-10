#pragma once

#include "../../Signal/Signal.hpp"
#include "../FilterUtility.hpp"
#include "../Windowing.hpp"

#include <type_traits>


namespace dspbb {

namespace impl {

	/// <summary> The default response for FIR filters, which is all-pass. </summary>
	struct {
		template <class T>
		auto operator()(T f) const { return T(1); }
	} inline constexpr DefaultResponse{};


	/// <summary> The default weights for how closely match desired and actual response for FIR filters, evenly weighted. </summary>
	struct {
		template <class T>
		auto operator()(T f) const { return T(1); }
	} inline constexpr DefaultWeight{};


	namespace windowed {

		template <template <typename, typename...> class Desc, class ParamType, class WindowType>
		struct SplitDescWindowed {
			ParamType cutoff = ParamType(0.5);
			WindowType window;

			/// <summary> Set the cutoff frequency of the filter. </summary>
			template <class NewParamType>
			[[nodiscard]] auto Cutoff(NewParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<NewParamType, WindowType>{ { std::move(cutoffNew), window } };
			}

			/// <summary> Set the window used for the window method when creating the filter as a generator function. </summary>
			template <windows_function_factory NewWindowType>
			[[nodiscard]] auto Window(NewWindowType windowNew) const {
				return Desc<ParamType, NewWindowType>{ { cutoff, std::move(windowNew) } };
			}

			/// <summary> Set the window used for the window method when creating the filter as a signal. </summary>
			template <signal_or_view NewWindowType>
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

			/// <summary> Set the band of the filter. </summary>
			template <class NewParamType>
			[[nodiscard]] auto Band(NewParamType lowerNew, NewParamType upperNew) const {
				impl::ThrowIfNotNormalized(lowerNew);
				impl::ThrowIfNotNormalized(upperNew);
				impl::ThrowIfNotSorted(lowerNew, upperNew);
				return Desc<NewParamType, WindowType>{ { std::move(lowerNew), std::move(upperNew), window } };
			}

			/// <summary> Set the window used for the window method when creating the filter as a generator function. </summary>
			template <windows_function_factory NewWindowType>
			[[nodiscard]] auto Window(NewWindowType windowNew) const {
				return Desc<ParamType, NewWindowType>{ { lower, upper, std::move(windowNew) } };
			}

			/// <summary> Set the window used for the window method when creating the filter as a signal. </summary>
			template <signal_or_view NewWindowType>
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

			template <std::invocable<float> NewResponseFunc>
			[[nodiscard]] auto Response(NewResponseFunc responseFuncNew) const {
				return ArbitraryDesc<NewResponseFunc, WindowType>{ std::move(responseFuncNew), window };
			}

			/// <summary> Set the window used for the window method when creating the filter as a generator function. </summary>
			template <windows_function_factory NewWindowType>
			[[nodiscard]] auto Window(NewWindowType windowNew) const {
				return ArbitraryDesc<ResponseFunc, NewWindowType>{ responseFunc, std::move(windowNew) };
			}

			/// <summary> Set the window used for the window method when creating the filter as a signal. </summary>
			template <signal_or_view NewWindowType>
			[[nodiscard]] auto Window(NewWindowType windowNew) const {
				return ArbitraryDesc<ResponseFunc, NewWindowType>{ responseFunc, std::move(windowNew) };
			}
		};

		template <class WindowType>
		struct HilbertDesc {
			WindowType window;

			/// <summary> Set the window used for the window method when creating the filter as a generator function. </summary>
			template <windows_function_factory NewWindowType>
			[[nodiscard]] auto Window(NewWindowType windowNew) const {
				return HilbertDesc<NewWindowType>{ std::move(windowNew) };
			}

			/// <summary> Set the window used for the window method when creating the filter as a signal. </summary>
			template <signal_or_view NewWindowType>
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

			/// <summary> Set the beginning and the end of the cutoff of the filter. </summary>
			/// <remarks> The response between the beginning and the end of the cutoff region
			///		will smoothly transition between stop and pass characteristics.
			///		A larger range will result in less overshoot with the same filter size. </remarks>
			[[nodiscard]] auto Cutoff(ParamType cutoffBeginNew, ParamType cutoffEndNew) const {
				impl::ThrowIfNotNormalized(cutoffBeginNew);
				impl::ThrowIfNotNormalized(cutoffEndNew);
				impl::ThrowIfNotSorted(cutoffBeginNew, cutoffEndNew);
				return Desc<ParamType>{ { cutoffBeginNew, cutoffEndNew, weightLow, weightTransition, weightHigh, grid } };
			}

			/// <summary> Set the relative importance of filter regions. </summary>
			///	<remarks> For a low-pass filter, a higher weight for the high region will
			///		result in better attenuation of high frequencies, but may add more
			///		ripples and less accurate amplification at the low frequencies. </remarks>
			[[nodiscard]] auto Weight(ParamType newLow, ParamType newTransition, ParamType newHigh) const {
				return Desc<ParamType>{ { cutoffBegin, cutoffEnd, newLow, newTransition, newHigh, grid } };
			}

			/// <summary> Set the number of points to which the response is discretized. </summary>
			/// <remarks> The least-squares method must discretize the desired continuous filter response.
			///		The finer the discretization, the less erratic the actual response, but
			///		the more expensive to compute the filter. </remarks>
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

			/// <summary> Set the beginning and the end of the transition regions of the filter. </summary>
			/// <remarks> Larger ranges will result in less overshoot with the same filter size. </remarks>
			[[nodiscard]] auto Band(ParamType lowerBeginNew, ParamType lowerEndNew, ParamType upperBeginNew, ParamType upperEndNew) const {
				impl::ThrowIfNotNormalized(lowerBeginNew);
				impl::ThrowIfNotNormalized(lowerEndNew);
				impl::ThrowIfNotNormalized(upperBeginNew);
				impl::ThrowIfNotNormalized(upperEndNew);
				impl::ThrowIfNotSorted(lowerBeginNew, lowerEndNew, upperBeginNew, upperEndNew);
				return Desc<ParamType>{ { lowerBeginNew, lowerEndNew, upperBeginNew, upperEndNew, weightLow, weightTransition1, weightMid, weightTransition2, weightHigh, grid } };
			}

			/// <summary> Set the relative importance of filter regions. </summary>
			///	<remarks> For a band-stop filter, a higher weight for the stop region will
			///		result in better attenuation of stop-frequencies, but may add more
			///		ripples and less accurate amplification at the pass regions. </remarks>
			[[nodiscard]] auto Weight(ParamType lowNew, ParamType transition1New, ParamType midNew, ParamType transition2New, ParamType highNew) const {
				return Desc<ParamType>{ { lowerBegin, lowerEnd, upperBegin, upperEnd, lowNew, transition1New, midNew, transition2New, highNew, grid } };
			}

			/// <summary> Set the number of points to which the response is discretized. </summary>
			/// <remarks> The least-squares method must discretize the desired continuous filter response.
			///		The finer the discretization, the less erratic the actual response, but
			///		the more expensive to compute the filter. </remarks>
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
			ParamType transitionWeight = ParamType(1.0);
			size_t grid = 0;

			/// <summary> Set the width of the transition region. </summary>
			/// <remarks> An FIR Hilbert filter's response tapers off near zero
			///		and the Nyquist frequency. The transition width refers to the width
			///		of the tapered ragion. </remarks>
			[[nodiscard]] auto TransitionWidth(ParamType newTransitionWidth) const {
				return HilbertDesc<ParamType>{ newTransitionWidth, transitionWeight, grid };
			}

			/// <summary> Set the weight of the transition region. </summary>
			/// <remarks> Setting this to zero might produce erratic behavior in the transition region.
			///		It's recommended to set it to something smaller than 1. </remarks>
			[[nodiscard]] auto TransitionWeight(ParamType newTransitionWeight) const {
				return HilbertDesc<ParamType>{ transitionWidth, newTransitionWeight, grid };
			}

			/// <summary> Set the number of points to which the response is discretized. </summary>
			/// <remarks> The least-squares method must discretize the desired continuous filter response.
			///		The finer the discretization, the less erratic the actual response, but
			///		the more expensive to compute the filter. </remarks>
			[[nodiscard]] auto Grid(size_t gridNew) const {
				return HilbertDesc<ParamType>{ transitionWidth, transitionWeight, gridNew };
			}
		};


		template <class ResponseFunc, class WeightFunc>
		struct ArbitraryDesc {
			ResponseFunc responseFunc{};
			WeightFunc weightFunc{};
			size_t grid = 0;

			/// <summary> Set the response of this filter. </summary>
			/// <param name="responseFuncNew"> Any real->real function. </param>
			template <std::invocable<float> NewResponseFunc>
			[[nodiscard]] auto Response(NewResponseFunc responseFuncNew) const {
				return ArbitraryDesc<NewResponseFunc, WeightFunc>{ std::move(responseFuncNew), weightFunc, grid };
			}

			/// <summary> Set the weights of this filter. </summary>
			/// <param name="weightFuncNew"> Any real->real function. </param>
			template <std::invocable<float> NewWeightFunc>
			[[nodiscard]] auto Weight(NewWeightFunc weightFuncNew) const {
				return ArbitraryDesc<ResponseFunc, NewWeightFunc>{ responseFunc, std::move(weightFuncNew), grid };
			}

			/// <summary> Set the number of points to which the response is discretized. </summary>
			/// <remarks> The least-squares method must discretize the desired continuous filter response.
			///		The finer the discretization, the less erratic the actual response, but
			///		the more expensive to compute the filter. </remarks>
			[[nodiscard]] auto Grid(size_t gridNew) const {
				return ArbitraryDesc<ResponseFunc, WeightFunc>{ responseFunc, weightFunc, gridNew };
			}
		};


		template <template <typename> class Desc>
		struct SplitDescLeastSquares<Desc, void> {
			/// <summary> Set the beginning and the end of the cutoff of the filter. </summary>
			/// <remarks> The response between the beginning and the end of the cutoff region
			///		will smoothly transition between stop and pass characteristics.
			///		A larger range will result in less overshoot with the same filter size. </remarks>
			template <class ParamType>
			[[nodiscard]] auto Cutoff(ParamType cutoffBeginNew, ParamType cutoffEndNew) const {
				return Desc<ParamType>{}.Cutoff(cutoffBeginNew, cutoffEndNew);
			}

			/// <summary> Set the relative importance of filter regions. </summary>
			///	<remarks> For a low-pass filter, a higher weight for the high region will
			///		result in better attenuation of high frequencies, but may add more
			///		ripples and less accurate amplification at the low frequencies. </remarks>
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
			/// <summary> Set the beginning and the end of the transition regions of the filter. </summary>
			/// <remarks> Larger ranges will result in less overshoot with the same filter size. </remarks>
			template <class ParamType>
			[[nodiscard]] auto Band(ParamType lowerBeginNew, ParamType lowerEndNew, ParamType upperBeginNew, ParamType upperEndNew) const {
				return Desc<ParamType>{}.Band(lowerBeginNew, lowerEndNew, upperBeginNew, upperEndNew);
			}

			/// <summary> Set the relative importance of filter regions. </summary>
			///	<remarks> For a band-stop filter, a higher weight for the stop region will
			///		result in better attenuation of stop-frequencies, but may add more
			///		ripples and less accurate amplification at the pass regions. </remarks>
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
			/// <summary> Set the width of the transition region. </summary>
			/// <remarks> An FIR Hilbert filter's response tapers off near zero
			///		and the Nyquist frequency. The transition width refers to the width
			///		of the tapered ragion. </remarks>
			template <class ParamType>
			[[nodiscard]] auto TransitionWidth(ParamType newTransitionWidth) const {
				return HilbertDesc<ParamType>{}.TransitionWidth(newTransitionWidth);
			}

			/// <summary> Set the weight of the transition region. </summary>
			/// <remarks> Setting this to zero might produce erratic behavior in the transition region.
			///		It's recommended to set it to something smaller than 1. </remarks>
			template <class ParamType>
			[[nodiscard]] auto TransitionWeight(ParamType newTransitionWeight) const {
				return HilbertDesc<ParamType>{}.TransitionWeight(newTransitionWeight);
			}
		};


	} // namespace least_squares

} // namespace impl

//------------------------------------------------------------------------------
// Factory functions
//------------------------------------------------------------------------------


/// <summary> FIR filter descriptions. </summary>
struct {
	/// <summary> Low-pass filter descriptions. </summary>
	struct {
		/// <summary> Description of a windowed low-pass filter. </summary>
		const impl::windowed::LowpassDesc<float, windows::Hamming> Windowed{};

		/// <summary> Description of a least-squares low-pass filter. </summary>
		const impl::least_squares::LowpassDesc<void> LeastSquares{};
	} const Lowpass{};

	/// <summary> High-pass filter descriptions. </summary>
	struct {
		/// <summary> Description of a windowed high-pass filter. </summary>
		const impl::windowed::HighpassDesc<float, windows::Hamming> Windowed{};

		/// <summary> Description of a least-squares high-pass filter. </summary>
		const impl::least_squares::HighpassDesc<void> LeastSquares{};
	} const Highpass{};
	struct {
		/// <summary> Description of a windowed band-pass filter. </summary>
		const impl::windowed::BandpassDesc<float, windows::Hamming> Windowed{};

		/// <summary> Description of a least-squares band-pass filter. </summary>
		const impl::least_squares::BandpassDesc<void> LeastSquares{};
	} const Bandpass{};
	struct {
		/// <summary> Description of a windowed band-stop filter. </summary>
		const impl::windowed::BandstopDesc<float, windows::Hamming> Windowed{};

		/// <summary> Description of a least-squares band-stop filter. </summary>
		const impl::least_squares::BandstopDesc<void> LeastSquares{};
	} const Bandstop{};
	struct {
		/// <summary> Description of a windowed Hilbert filter. </summary>
		const impl::windowed::HilbertDesc<windows::Hamming> Windowed{};

		/// <summary> Description of a least-squares Hilbert filter. </summary>
		const impl::least_squares::HilbertDesc<void> LeastSquares{};
	} const Hilbert{};
	struct {
		/// <summary> Description of a windowed arbitrary response filter. </summary>
		const impl::windowed::ArbitraryDesc<decltype(impl::DefaultResponse), windows::Hamming> Windowed{};

		/// <summary> Description of a least-squares arbitrary response filter. </summary>
		const impl::least_squares::ArbitraryDesc<decltype(impl::DefaultResponse), decltype(impl::DefaultWeight)> LeastSquares{};
	} const Arbitrary{};
} inline const Fir;


} // namespace dspbb
