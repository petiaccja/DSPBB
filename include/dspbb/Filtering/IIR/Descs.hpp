#pragma once

#include "../FilterUtility.hpp"

#include <type_traits>
#include <utility>


namespace dspbb {


namespace impl {

	namespace butterworth {

		template <template <typename> class Desc, class ParamType>
		struct SplitDescBase {
			ParamType cutoff = ParamType(0.5);

			/// <summary> Set the cutoff frequency of the filter, in normalized frequency. </summary>
			template <class NewParamType>
			[[nodiscard]] auto Cutoff(NewParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<NewParamType>{ { cutoffNew } };
			}
		};

		template <template <typename> class Desc, class ParamType>
		struct BandDescBase {
			ParamType lower = ParamType(0.25);
			ParamType upper = ParamType(0.75);

			/// <summary> Set the beginning and end frequency of the filter's band, in normalized frequencies. </summary>
			template <class NewParamType>
			[[nodiscard]] auto Band(NewParamType lowerNew, NewParamType upperNew) const {
				impl::ThrowIfNotNormalized(lowerNew);
				impl::ThrowIfNotNormalized(upperNew);
				impl::ThrowIfNotSorted(lowerNew, upperNew);
				return Desc<NewParamType>{ { lowerNew, upperNew } };
			}
		};

		template <class T>
		struct LowpassDesc : SplitDescBase<LowpassDesc, T> {};

		template <class T>
		struct HighpassDesc : SplitDescBase<HighpassDesc, T> {};

		template <class T>
		struct BandpassDesc : BandDescBase<BandpassDesc, T> {};

		template <class T>
		struct BandstopDesc : BandDescBase<BandstopDesc, T> {};

	} // namespace butterworth


	namespace chebyshev1 {

		template <template <typename> class Desc, class ParamType>
		struct SplitDescBase {
			ParamType cutoff = ParamType(0.5);
			ParamType passbandRipple = ParamType(0.1);

			/// <summary> Set the cutoff frequency of the filter, in normalized frequency. </summary>
			[[nodiscard]] auto Cutoff(ParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<ParamType>{ { cutoffNew, passbandRipple } };
			}

			/// <summary> Set the magnitude of the passband ripple. </summary>
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { cutoff, rippleNew } };
			}
		};

		template <template <typename> class Desc, class ParamType>
		struct BandDescBase {
			ParamType lower = ParamType(0.25);
			ParamType upper = ParamType(0.75);
			ParamType passbandRipple = ParamType(0.1);

			/// <summary> Set the beginning and end frequency of the filter's band, in normalized frequencies. </summary>
			[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) const {
				impl::ThrowIfNotNormalized(lowerNew);
				impl::ThrowIfNotNormalized(upperNew);
				impl::ThrowIfNotSorted(lowerNew, upperNew);
				return Desc<ParamType>{ { lowerNew, upperNew, passbandRipple } };
			}

			/// <summary> Set the magnitude of the passband ripple. </summary>
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { lower, upper, rippleNew } };
			}
		};

		template <template <typename> class Desc>
		struct SplitDescBase<Desc, void> {
			/// <summary> Set the cutoff frequency of the filter, in normalized frequency. </summary>
			template <class ParamType>
			[[nodiscard]] auto Cutoff(ParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<ParamType>{}.Cutoff(cutoffNew);
			}

			/// <summary> Set the magnitude of the passband ripple. </summary>
			template <class ParamType>
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{}.PassbandRipple(rippleNew);
			}
		};

		template <template <typename> class Desc>
		struct BandDescBase<Desc, void> {
			/// <summary> Set the beginning and end frequency of the filter's band, in normalized frequencies. </summary>
			template <class ParamType>
			[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) const {
				impl::ThrowIfNotNormalized(lowerNew);
				impl::ThrowIfNotNormalized(upperNew);
				impl::ThrowIfNotSorted(lowerNew, upperNew);
				return Desc<ParamType>{}.Band(lowerNew, upperNew);
			}

			/// <summary> Set the magnitude of the passband ripple. </summary>
			template <class ParamType>
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{}.PassbandRipple(rippleNew);
			}
		};

		template <class T>
		struct LowpassDesc : SplitDescBase<LowpassDesc, T> {};

		template <>
		struct LowpassDesc<void> : SplitDescBase<LowpassDesc, void> {};

		template <class T>
		struct HighpassDesc : SplitDescBase<HighpassDesc, T> {};

		template <>
		struct HighpassDesc<void> : SplitDescBase<HighpassDesc, void> {};

		template <class T>
		struct BandpassDesc : BandDescBase<BandpassDesc, T> {};

		template <>
		struct BandpassDesc<void> : BandDescBase<BandpassDesc, void> {};

		template <class T>
		struct BandstopDesc : BandDescBase<BandstopDesc, T> {};

		template <>
		struct BandstopDesc<void> : BandDescBase<BandstopDesc, void> {};

	} // namespace chebyshev1


	namespace chebyshev2 {

		template <template <typename> class Desc, class ParamType>
		struct SplitDescBase {
			ParamType cutoff = ParamType(0.5);
			ParamType stopbandRipple = ParamType(0.1);

			/// <summary> Set the cutoff frequency of the filter, in normalized frequency. </summary>
			[[nodiscard]] auto Cutoff(ParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<ParamType>{ { cutoffNew, stopbandRipple } };
			}

			/// <summary> Set the magnitude of the stopband ripple. </summary>
			[[nodiscard]] auto StopbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { cutoff, rippleNew } };
			}
		};

		template <template <typename> class Desc, class ParamType>
		struct BandDescBase {
			ParamType lower = ParamType(0.25);
			ParamType upper = ParamType(0.75);
			ParamType stopbandRipple = ParamType(0.1);

			/// <summary> Set the beginning and end frequency of the filter's band, in normalized frequencies. </summary>
			[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) const {
				impl::ThrowIfNotNormalized(lowerNew);
				impl::ThrowIfNotNormalized(upperNew);
				impl::ThrowIfNotSorted(lowerNew, upperNew);
				return Desc<ParamType>{ { lowerNew, upperNew, stopbandRipple } };
			}

			/// <summary> Set the magnitude of the stopband ripple. </summary>
			[[nodiscard]] auto StopbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { lower, upper, rippleNew } };
			}
		};

		template <template <typename> class Desc>
		struct SplitDescBase<Desc, void> {
			/// <summary> Set the cutoff frequency of the filter, in normalized frequency. </summary>
			template <class ParamType>
			[[nodiscard]] auto Cutoff(ParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<ParamType>{}.Cutoff(cutoffNew);
			}

			/// <summary> Set the magnitude of the stopband ripple. </summary>
			template <class ParamType>
			[[nodiscard]] auto StopbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{}.StopbandRipple(rippleNew);
			}
		};

		template <template <typename> class Desc>
		struct BandDescBase<Desc, void> {
			/// <summary> Set the beginning and end frequency of the filter's band, in normalized frequencies. </summary>
			template <class ParamType>
			[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) const {
				impl::ThrowIfNotNormalized(lowerNew);
				impl::ThrowIfNotNormalized(upperNew);
				impl::ThrowIfNotSorted(lowerNew, upperNew);
				return Desc<ParamType>{}.Band(lowerNew, upperNew);
			}

			/// <summary> Set the magnitude of the stopband ripple. </summary>
			template <class ParamType>
			[[nodiscard]] auto StopbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{}.StopbandRipple(rippleNew);
			}
		};

		template <class T>
		struct LowpassDesc : SplitDescBase<LowpassDesc, T> {};

		template <>
		struct LowpassDesc<void> : SplitDescBase<LowpassDesc, void> {};

		template <class T>
		struct HighpassDesc : SplitDescBase<HighpassDesc, T> {};

		template <>
		struct HighpassDesc<void> : SplitDescBase<HighpassDesc, void> {};

		template <class T>
		struct BandpassDesc : BandDescBase<BandpassDesc, T> {};

		template <>
		struct BandpassDesc<void> : BandDescBase<BandpassDesc, void> {};

		template <class T>
		struct BandstopDesc : BandDescBase<BandstopDesc, T> {};

		template <>
		struct BandstopDesc<void> : BandDescBase<BandstopDesc, void> {};
	} // namespace chebyshev2


	namespace elliptic {

		template <template <typename, typename...> class Desc, class ParamType>
		struct SplitDescBase {
			ParamType cutoff = ParamType(0.5);
			ParamType passbandRipple = ParamType(0.1);
			ParamType stopbandRipple = ParamType(0.1);

			/// <summary> Set the cutoff frequency of the filter, in normalized frequency. </summary>
			[[nodiscard]] auto Cutoff(ParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<ParamType>{ { cutoffNew, passbandRipple, stopbandRipple } };
			}

			/// <summary> Set the magnitude of the passband ripple. </summary>
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { cutoff, rippleNew, stopbandRipple } };
			}

			/// <summary> Set the magnitude of the stopband ripple. </summary>
			[[nodiscard]] auto StopbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { cutoff, passbandRipple, rippleNew } };
			}
		};

		template <template <typename, typename...> class Desc, class ParamType>
		struct BandDescBase {
			ParamType lower = ParamType(0.25);
			ParamType upper = ParamType(0.75);
			ParamType passbandRipple = ParamType(0.1);
			ParamType stopbandRipple = ParamType(0.1);

			/// <summary> Set the beginning and end frequency of the filter's band, in normalized frequencies. </summary>
			[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) const {
				impl::ThrowIfNotNormalized(lowerNew);
				impl::ThrowIfNotNormalized(upperNew);
				impl::ThrowIfNotSorted(lowerNew, upperNew);
				return Desc<ParamType>{ { lowerNew, upperNew, passbandRipple, stopbandRipple } };
			}

			/// <summary> Set the magnitude of the passband ripple. </summary>
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { lower, upper, rippleNew, stopbandRipple } };
			}

			/// <summary> Set the magnitude of the stopband ripple. </summary>
			[[nodiscard]] auto StopbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { lower, upper, passbandRipple, rippleNew } };
			}
		};

		template <template <typename, typename...> class Desc>
		struct SplitDescBase<Desc, void> {
			/// <summary> Set the cutoff frequency of the filter, in normalized frequency. </summary>
			template <class ParamType>
			[[nodiscard]] auto Cutoff(ParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<ParamType>{}.Cutoff(cutoffNew);
			}

			/// <summary> Set the magnitude of the passband ripple. </summary>
			template <class ParamType>
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{}.PassbandRipple(rippleNew);
			}

			/// <summary> Set the magnitude of the stopband ripple. </summary>
			template <class ParamType>
			[[nodiscard]] auto StopbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{}.StopbandRipple(rippleNew);
			}
		};

		template <template <typename, typename...> class Desc>
		struct BandDescBase<Desc, void> {
			/// <summary> Set the beginning and end frequency of the filter's band, in normalized frequencies. </summary>
			template <class ParamType>
			[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) const {
				impl::ThrowIfNotNormalized(lowerNew);
				impl::ThrowIfNotNormalized(upperNew);
				impl::ThrowIfNotSorted(lowerNew, upperNew);
				return Desc<ParamType>{}.Band(lowerNew, upperNew);
			}

			/// <summary> Set the magnitude of the passband ripple. </summary>
			template <class ParamType>
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{}.PassbandRipple(rippleNew);
			}

			/// <summary> Set the magnitude of the stopband ripple. </summary>
			template <class ParamType>
			[[nodiscard]] auto StopbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{}.StopbandRipple(rippleNew);
			}
		};

		template <class T>
		struct LowpassDesc : SplitDescBase<LowpassDesc, T> {};

		template <>
		struct LowpassDesc<void> : SplitDescBase<LowpassDesc, void> {};

		template <class T>
		struct HighpassDesc : SplitDescBase<HighpassDesc, T> {};

		template <>
		struct HighpassDesc<void> : SplitDescBase<HighpassDesc, void> {};

		template <class T>
		struct BandpassDesc : BandDescBase<BandpassDesc, T> {};

		template <>
		struct BandpassDesc<void> : BandDescBase<BandpassDesc, void> {};

		template <class T>
		struct BandstopDesc : BandDescBase<BandstopDesc, T> {};

		template <>
		struct BandstopDesc<void> : BandDescBase<BandstopDesc, void> {};

	} // namespace elliptic

} // namespace impl


/// <summary> IIR filter descriptions. </summary>
struct {
	/// <summary> Low-pass filter descriptions. </summary>
	struct {
		/// <summary> Description of a Butterworth low-pass filter. </summary>
		const impl::butterworth::LowpassDesc<float> Butterworth{};
		/// <summary> Description of a Chebyshev type 1 low-pass filter. </summary>
		const impl::chebyshev1::LowpassDesc<void> Chebyshev1{};
		/// <summary> Description of a Chebyshev type 2 low-pass filter. </summary>
		const impl::chebyshev2::LowpassDesc<void> Chebyshev2{};
		/// <summary> Description of a elliptic low-pass filter. </summary>
		const impl::elliptic::LowpassDesc<void> Elliptic{};
	} Lowpass{};
	/// <summary> High-pass filter descriptions. </summary>
	struct {
		/// <summary> Description of a Butterworth high-pass filter. </summary>
		const impl::butterworth::HighpassDesc<float> Butterworth{};
		/// <summary> Description of a Chebyshev type 1 high-pass filter. </summary>
		const impl::chebyshev1::HighpassDesc<void> Chebyshev1{};
		/// <summary> Description of a Chebyshev type 2 high-pass filter. </summary>
		const impl::chebyshev2::HighpassDesc<void> Chebyshev2{};
		/// <summary> Description of a elliptic high-pass filter. </summary>
		const impl::elliptic::HighpassDesc<void> Elliptic{};
	} Highpass{};
	/// <summary> Band-pass filter descriptions. </summary>
	struct {
		/// <summary> Description of a Butterworth band-pass filter. </summary>
		const impl::butterworth::BandpassDesc<float> Butterworth{};
		/// <summary> Description of a Chebyshev type 1 band-pass filter. </summary>
		const impl::chebyshev1::BandpassDesc<void> Chebyshev1{};
		/// <summary> Description of a Chebyshev type 2 band-pass filter. </summary>
		const impl::chebyshev2::BandpassDesc<void> Chebyshev2{};
		/// <summary> Description of a elliptic band-pass filter. </summary>
		const impl::elliptic::BandpassDesc<void> Elliptic{};
	} Bandpass{};
	/// <summary> Band-stop filter descriptions. </summary>
	struct {
		/// <summary> Description of a Butterworth band-stopelliptic filter. </summary>
		const impl::butterworth::BandstopDesc<float> Butterworth{};
		/// <summary> Description of a Chebyshev type 1 band-stop filter. </summary>
		const impl::chebyshev1::BandstopDesc<void> Chebyshev1{};
		/// <summary> Description of a Chebyshev type 2 band-stop filter. </summary>
		const impl::chebyshev2::BandstopDesc<void> Chebyshev2{};
		/// <summary> Description of a elliptic band-stop filter. </summary>
		const impl::elliptic::BandstopDesc<void> Elliptic{};
	} Bandstop{};
} const Iir{};


} // namespace dspbb