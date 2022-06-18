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

			[[nodiscard]] auto Cutoff(ParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<ParamType>{ { cutoffNew, passbandRipple } };
			}
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { cutoff, rippleNew } };
			}
		};

		template <template <typename> class Desc, class ParamType>
		struct BandDescBase {
			ParamType lower = ParamType(0.25);
			ParamType upper = ParamType(0.75);
			ParamType passbandRipple = ParamType(0.1);

			[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) const {
				impl::ThrowIfNotNormalized(lowerNew);
				impl::ThrowIfNotNormalized(upperNew);
				impl::ThrowIfNotSorted(lowerNew, upperNew);
				return Desc<ParamType>{ { lowerNew, upperNew, passbandRipple } };
			}
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { lower, upper, rippleNew } };
			}
		};

		template <template <typename> class Desc>
		struct SplitDescBase<Desc, void> {
			template <class ParamType>
			[[nodiscard]] auto Cutoff(ParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<ParamType>{}.Cutoff(cutoffNew);
			}
			template <class ParamType>
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{}.PassbandRipple(rippleNew);
			}
		};

		template <template <typename> class Desc>
		struct BandDescBase<Desc, void> {
			template <class ParamType>
			[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) const {
				impl::ThrowIfNotNormalized(lowerNew);
				impl::ThrowIfNotNormalized(upperNew);
				impl::ThrowIfNotSorted(lowerNew, upperNew);
				return Desc<ParamType>{}.Band(lowerNew, upperNew);
			}
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

			[[nodiscard]] auto Cutoff(ParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<ParamType>{ { cutoffNew, stopbandRipple } };
			}
			[[nodiscard]] auto StopbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { cutoff, rippleNew } };
			}
		};

		template <template <typename> class Desc, class ParamType>
		struct BandDescBase {
			ParamType lower = ParamType(0.25);
			ParamType upper = ParamType(0.75);
			ParamType stopbandRipple = ParamType(0.1);

			[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) const {
				impl::ThrowIfNotNormalized(lowerNew);
				impl::ThrowIfNotNormalized(upperNew);
				impl::ThrowIfNotSorted(lowerNew, upperNew);
				return Desc<ParamType>{ { lowerNew, upperNew, stopbandRipple } };
			}
			[[nodiscard]] auto StopbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { lower, upper, rippleNew } };
			}
		};

		template <template <typename> class Desc>
		struct SplitDescBase<Desc, void> {
			template <class ParamType>
			[[nodiscard]] auto Cutoff(ParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<ParamType>{}.Cutoff(cutoffNew);
			}
			template <class ParamType>
			[[nodiscard]] auto StopbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{}.StopbandRipple(rippleNew);
			}
		};

		template <template <typename> class Desc>
		struct BandDescBase<Desc, void> {
			template <class ParamType>
			[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) const {
				impl::ThrowIfNotNormalized(lowerNew);
				impl::ThrowIfNotNormalized(upperNew);
				impl::ThrowIfNotSorted(lowerNew, upperNew);
				return Desc<ParamType>{}.Band(lowerNew, upperNew);
			}
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

			[[nodiscard]] auto Cutoff(ParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<ParamType>{ { cutoffNew, passbandRipple, stopbandRipple } };
			}
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { cutoff, rippleNew, stopbandRipple } };
			}
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

			[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) const {
				impl::ThrowIfNotNormalized(lowerNew);
				impl::ThrowIfNotNormalized(upperNew);
				impl::ThrowIfNotSorted(lowerNew, upperNew);
				return Desc<ParamType>{ { lowerNew, upperNew, passbandRipple, stopbandRipple } };
			}
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { lower, upper, rippleNew, stopbandRipple } };
			}
			[[nodiscard]] auto StopbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{ { lower, upper, passbandRipple, rippleNew } };
			}
		};

		template <template <typename, typename...> class Desc>
		struct SplitDescBase<Desc, void> {
			template <class ParamType>
			[[nodiscard]] auto Cutoff(ParamType cutoffNew) const {
				impl::ThrowIfNotNormalized(cutoffNew);
				return Desc<ParamType>{}.Cutoff(cutoffNew);
			}
			template <class ParamType>
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{}.PassbandRipple(rippleNew);
			}
			template <class ParamType>
			[[nodiscard]] auto StopbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{}.StopbandRipple(rippleNew);
			}
		};

		template <template <typename, typename...> class Desc>
		struct BandDescBase<Desc, void> {
			template <class ParamType>
			[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) const {
				impl::ThrowIfNotNormalized(lowerNew);
				impl::ThrowIfNotNormalized(upperNew);
				impl::ThrowIfNotSorted(lowerNew, upperNew);
				return Desc<ParamType>{}.Band(lowerNew, upperNew);
			}
			template <class ParamType>
			[[nodiscard]] auto PassbandRipple(ParamType rippleNew) const {
				return Desc<ParamType>{}.PassbandRipple(rippleNew);
			}
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


struct {
	struct {
		const impl::butterworth::LowpassDesc<float> Butterworth{};
		const impl::chebyshev1::LowpassDesc<void> Chebyshev1{};
		const impl::chebyshev2::LowpassDesc<void> Chebyshev2{};
		const impl::elliptic::LowpassDesc<void> Elliptic{};
	} Lowpass{};
	struct {
		const impl::butterworth::HighpassDesc<float> Butterworth{};
		const impl::chebyshev1::HighpassDesc<void> Chebyshev1{};
		const impl::chebyshev2::HighpassDesc<void> Chebyshev2{};
		const impl::elliptic::HighpassDesc<void> Elliptic{};
	} Highpass{};
	struct {
		const impl::butterworth::BandpassDesc<float> Butterworth{};
		const impl::chebyshev1::BandpassDesc<void> Chebyshev1{};
		const impl::chebyshev2::BandpassDesc<void> Chebyshev2{};
		const impl::elliptic::BandpassDesc<void> Elliptic{};
	} Bandpass{};
	struct {
		const impl::butterworth::BandstopDesc<float> Butterworth{};
		const impl::chebyshev1::BandstopDesc<void> Chebyshev1{};
		const impl::chebyshev2::BandstopDesc<void> Chebyshev2{};
		const impl::elliptic::BandstopDesc<void> Elliptic{};
	} Bandstop{};
} const Iir{};


} // namespace dspbb