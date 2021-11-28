#pragma once

#include <type_traits>
#include <utility>


namespace dspbb {


namespace impl {
	struct IirMethod {};
	struct IirMethodButterworth : IirMethod {};
	struct IirMethodChebyshev1 : IirMethod {};
	struct IirMethodChebyshev2 : IirMethod {};
	struct IirMethodElliptic : IirMethod {};


	template <class Method, class... Params>
	struct LowpassDesc;

	template <class Method, class... Params>
	struct HighpassDesc;

	template <class Method, class... Params>
	struct BandpassDesc;

	template <class Method, class... Params>
	struct BandstopDesc;

	//------------------------------------------------------------------------------
	// Butterworth descs
	//------------------------------------------------------------------------------

	template <template <typename, typename...> class Desc, class ParamType>
	struct SplitDescButterworth {
		ParamType cutoff = ParamType(0.5);

		template <class NewParamType>
		[[nodiscard]] auto Cutoff(NewParamType cutoffNew) {
			return Desc<IirMethodButterworth, NewParamType>{ { std::move(cutoffNew) } };
		}
	};

	template <template <typename, typename...> class Desc, class ParamType>
	struct BandDescButterworth {
		ParamType lower = ParamType(0.25);
		ParamType upper = ParamType(0.75);

		template <class NewParamType>
		[[nodiscard]] auto Band(NewParamType lowerNew, NewParamType upperNew) {
			return Desc<IirMethodButterworth, NewParamType>{ { std::move(lowerNew), std::move(upperNew) } };
		}
	};

	template <class T>
	struct LowpassDesc<IirMethodButterworth, T>
		: SplitDescButterworth<LowpassDesc, T> {};

	template <>
	struct LowpassDesc<IirMethodButterworth>
		: SplitDescButterworth<LowpassDesc, float> {};

	template <class T>
	struct HighpassDesc<IirMethodButterworth, T>
		: SplitDescButterworth<HighpassDesc, T> {};

	template <>
	struct HighpassDesc<IirMethodButterworth>
		: SplitDescButterworth<HighpassDesc, float> {};

	template <class T>
	struct BandpassDesc<IirMethodButterworth, T>
		: BandDescButterworth<BandpassDesc, T> {};

	template <>
	struct BandpassDesc<IirMethodButterworth>
		: BandDescButterworth<BandpassDesc, float> {};

	template <class T>
	struct BandstopDesc<IirMethodButterworth, T>
		: BandDescButterworth<BandstopDesc, T> {};

	template <>
	struct BandstopDesc<IirMethodButterworth>
		: BandDescButterworth<BandstopDesc, float> {};


	//------------------------------------------------------------------------------
	// Chebyshev 1 descs
	//------------------------------------------------------------------------------

	template <template <typename, typename...> class Desc, class ParamType>
	struct SplitDescChebyshev1 {
		ParamType cutoff = ParamType(0.5);
		ParamType passbandRipple = ParamType(0.1);

		[[nodiscard]] auto Cutoff(ParamType cutoffNew) {
			return Desc<IirMethodChebyshev1, ParamType>{ { cutoffNew, passbandRipple } };
		}
		[[nodiscard]] auto PassbandRipple(ParamType rippleNew) {
			return Desc<IirMethodChebyshev1, ParamType>{ { cutoff, rippleNew } };
		}
	};

	template <template <typename, typename...> class Desc, class ParamType>
	struct BandDescChebyshev1 {
		ParamType lower = ParamType(0.25);
		ParamType upper = ParamType(0.75);
		ParamType passbandRipple = ParamType(0.1);

		[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) {
			return Desc<IirMethodChebyshev1, ParamType>{ { lowerNew, upperNew, passbandRipple } };
		}
		[[nodiscard]] auto PassbandRipple(ParamType rippleNew) {
			return Desc<IirMethodChebyshev1, ParamType>{ { lower, upper, rippleNew } };
		}
	};

	template <template <typename, typename...> class Desc>
	struct SplitDescChebyshev1<Desc, void> {
		template <class ParamType>
		[[nodiscard]] auto Cutoff(ParamType cutoffNew) {
			return Desc<IirMethodChebyshev1, ParamType>{}.Cutoff(cutoffNew);
		}
		template <class ParamType>
		[[nodiscard]] auto PassbandRipple(ParamType rippleNew) {
			return Desc<IirMethodChebyshev1, ParamType>{}.PassbandRipple(rippleNew);
		}
	};

	template <template <typename, typename...> class Desc>
	struct BandDescChebyshev1<Desc, void> {
		template <class ParamType>
		[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) {
			return Desc<IirMethodChebyshev1, ParamType>{}.Band(lowerNew, upperNew);
		}
		template <class ParamType>
		[[nodiscard]] auto PassbandRipple(ParamType rippleNew) {
			return Desc<IirMethodChebyshev1, ParamType>{}.PassbandRipple(rippleNew);
		}
	};

	template <class T>
	struct LowpassDesc<IirMethodChebyshev1, T>
		: SplitDescChebyshev1<LowpassDesc, T> {};

	template <>
	struct LowpassDesc<IirMethodChebyshev1>
		: SplitDescChebyshev1<LowpassDesc, void> {};

	template <class T>
	struct HighpassDesc<IirMethodChebyshev1, T>
		: SplitDescChebyshev1<HighpassDesc, T> {};

	template <>
	struct HighpassDesc<IirMethodChebyshev1>
		: SplitDescChebyshev1<HighpassDesc, void> {};

	template <class T>
	struct BandpassDesc<IirMethodChebyshev1, T>
		: BandDescChebyshev1<BandpassDesc, T> {};

	template <>
	struct BandpassDesc<IirMethodChebyshev1>
		: BandDescChebyshev1<BandpassDesc, void> {};

	template <class T>
	struct BandstopDesc<IirMethodChebyshev1, T>
		: BandDescChebyshev1<BandstopDesc, T> {};

	template <>
	struct BandstopDesc<IirMethodChebyshev1>
		: BandDescChebyshev1<BandstopDesc, void> {};

	//------------------------------------------------------------------------------
	// Chebyshev 2 descs
	//------------------------------------------------------------------------------

	template <template <typename, typename...> class Desc, class ParamType>
	struct SplitDescChebyshev2 {
		ParamType cutoff = ParamType(0.5);
		ParamType stopbandRipple = ParamType(0.1);

		[[nodiscard]] auto Cutoff(ParamType cutoffNew) {
			return Desc<IirMethodChebyshev2, ParamType>{ { cutoffNew, stopbandRipple } };
		}
		[[nodiscard]] auto StopbandRipple(ParamType rippleNew) {
			return Desc<IirMethodChebyshev2, ParamType>{ { cutoff, rippleNew } };
		}
	};

	template <template <typename, typename...> class Desc, class ParamType>
	struct BandDescChebyshev2 {
		ParamType lower = ParamType(0.25);
		ParamType upper = ParamType(0.75);
		ParamType stopbandRipple = ParamType(0.1);

		[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) {
			return Desc<IirMethodChebyshev2, ParamType>{ { lowerNew, upperNew, stopbandRipple } };
		}
		[[nodiscard]] auto StopbandRipple(ParamType rippleNew) {
			return Desc<IirMethodChebyshev2, ParamType>{ { lower, upper, rippleNew } };
		}
	};

	template <template <typename, typename...> class Desc>
	struct SplitDescChebyshev2<Desc, void> {
		template <class ParamType>
		[[nodiscard]] auto Cutoff(ParamType cutoffNew) {
			return Desc<IirMethodChebyshev2, ParamType>{}.Cutoff(cutoffNew);
		}
		template <class ParamType>
		[[nodiscard]] auto StopbandRipple(ParamType rippleNew) {
			return Desc<IirMethodChebyshev2, ParamType>{}.StopbandRipple(rippleNew);
		}
	};

	template <template <typename, typename...> class Desc>
	struct BandDescChebyshev2<Desc, void> {
		template <class ParamType>
		[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) {
			return Desc<IirMethodChebyshev2, ParamType>{}.Band(lowerNew, upperNew);
		}
		template <class ParamType>
		[[nodiscard]] auto StopbandRipple(ParamType rippleNew) {
			return Desc<IirMethodChebyshev2, ParamType>{}.StopbandRipple(rippleNew);
		}
	};

	template <class T>
	struct LowpassDesc<IirMethodChebyshev2, T>
		: SplitDescChebyshev2<LowpassDesc, T> {};

	template <>
	struct LowpassDesc<IirMethodChebyshev2>
		: SplitDescChebyshev2<LowpassDesc, void> {};

	template <class T>
	struct HighpassDesc<IirMethodChebyshev2, T>
		: SplitDescChebyshev2<HighpassDesc, T> {};

	template <>
	struct HighpassDesc<IirMethodChebyshev2>
		: SplitDescChebyshev2<HighpassDesc, void> {};

	template <class T>
	struct BandpassDesc<IirMethodChebyshev2, T>
		: BandDescChebyshev2<BandpassDesc, T> {};

	template <>
	struct BandpassDesc<IirMethodChebyshev2>
		: BandDescChebyshev2<BandpassDesc, void> {};

	template <class T>
	struct BandstopDesc<IirMethodChebyshev2, T>
		: BandDescChebyshev2<BandstopDesc, T> {};

	template <>
	struct BandstopDesc<IirMethodChebyshev2>
		: BandDescChebyshev2<BandstopDesc, void> {};

	//------------------------------------------------------------------------------
	// Elliptic descs
	//------------------------------------------------------------------------------


	template <template <typename, typename...> class Desc, class ParamType>
	struct SplitDescElliptic {
		ParamType cutoff = ParamType(0.5);
		ParamType passbandRipple = ParamType(0.1);
		ParamType stopbandRipple = ParamType(0.1);

		[[nodiscard]] auto Cutoff(ParamType cutoffNew) {
			return Desc<IirMethodElliptic, ParamType>{ { cutoffNew, passbandRipple, stopbandRipple } };
		}
		[[nodiscard]] auto PassbandRipple(ParamType rippleNew) {
			return Desc<IirMethodElliptic, ParamType>{ { cutoff, rippleNew, stopbandRipple } };
		}
		[[nodiscard]] auto StopbandRipple(ParamType rippleNew) {
			return Desc<IirMethodElliptic, ParamType>{ { cutoff, passbandRipple, rippleNew } };
		}
	};

	template <template <typename, typename...> class Desc, class ParamType>
	struct BandDescElliptic {
		ParamType lower = ParamType(0.25);
		ParamType upper = ParamType(0.75);
		ParamType passbandRipple = ParamType(0.1);
		ParamType stopbandRipple = ParamType(0.1);

		[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) {
			return Desc<IirMethodElliptic, ParamType>{ { lowerNew, upperNew, passbandRipple, stopbandRipple } };
		}
		[[nodiscard]] auto PassbandRipple(ParamType rippleNew) {
			return Desc<IirMethodElliptic, ParamType>{ { lower, upper, rippleNew, stopbandRipple } };
		}
		[[nodiscard]] auto StopbandRipple(ParamType rippleNew) {
			return Desc<IirMethodElliptic, ParamType>{ { lower, upper, passbandRipple, rippleNew } };
		}
	};

	template <template <typename, typename...> class Desc>
	struct SplitDescElliptic<Desc, void> {
		template <class ParamType>
		[[nodiscard]] auto Cutoff(ParamType cutoffNew) {
			return Desc<IirMethodElliptic, ParamType>{}.Cutoff(cutoffNew);
		}
		template <class ParamType>
		[[nodiscard]] auto PassbandRipple(ParamType rippleNew) {
			return Desc<IirMethodElliptic, ParamType>{}.PassbandRipple(rippleNew);
		}
		template <class ParamType>
		[[nodiscard]] auto StopbandRipple(ParamType rippleNew) {
			return Desc<IirMethodElliptic, ParamType>{}.StopbandRipple(rippleNew);
		}
	};

	template <template <typename, typename...> class Desc>
	struct BandDescElliptic<Desc, void> {
		template <class ParamType>
		[[nodiscard]] auto Band(ParamType lowerNew, ParamType upperNew) {
			return Desc<IirMethodElliptic, ParamType>{}.Band(lowerNew, upperNew);
		}
		template <class ParamType>
		[[nodiscard]] auto PassbandRipple(ParamType rippleNew) {
			return Desc<IirMethodElliptic, ParamType>{}.PassbandRipple(rippleNew);
		}
		template <class ParamType>
		[[nodiscard]] auto StopbandRipple(ParamType rippleNew) {
			return Desc<IirMethodElliptic, ParamType>{}.StopbandRipple(rippleNew);
		}
	};

	template <class T>
	struct LowpassDesc<IirMethodElliptic, T>
		: SplitDescElliptic<LowpassDesc, T> {};

	template <>
	struct LowpassDesc<IirMethodElliptic>
		: SplitDescElliptic<LowpassDesc, void> {};

	template <class T>
	struct HighpassDesc<IirMethodElliptic, T>
		: SplitDescElliptic<HighpassDesc, T> {};

	template <>
	struct HighpassDesc<IirMethodElliptic>
		: SplitDescElliptic<HighpassDesc, void> {};

	template <class T>
	struct BandpassDesc<IirMethodElliptic, T>
		: BandDescElliptic<BandpassDesc, T> {};

	template <>
	struct BandpassDesc<IirMethodElliptic>
		: BandDescElliptic<BandpassDesc, void> {};

	template <class T>
	struct BandstopDesc<IirMethodElliptic, T>
		: BandDescElliptic<BandstopDesc, T> {};

	template <>
	struct BandstopDesc<IirMethodElliptic>
		: BandDescElliptic<BandstopDesc, void> {};

} // namespace impl


constexpr impl::IirMethodButterworth BUTTERWORTH;
constexpr impl::IirMethodChebyshev1 CHEBYSHEV1;
constexpr impl::IirMethodChebyshev2 CHEBYSHEV2;
constexpr impl::IirMethodElliptic ELLIPTIC;


template <class Method, std::enable_if_t<std::is_base_of_v<impl::IirMethod, Method>, int> = 0>
auto Lowpass(Method) {
	return impl::LowpassDesc<Method>{};
}

template <class Method, std::enable_if_t<std::is_base_of_v<impl::IirMethod, Method>, int> = 0>
auto Highpass(Method) {
	return impl::HighpassDesc<Method>{};
}

template <class Method, std::enable_if_t<std::is_base_of_v<impl::IirMethod, Method>, int> = 0>
auto Bandpass(Method) {
	return impl::BandpassDesc<Method>{};
}

template <class Method, std::enable_if_t<std::is_base_of_v<impl::IirMethod, Method>, int> = 0>
auto Bandstop(Method) {
	return impl::BandstopDesc<Method>{};
}


} // namespace dspbb