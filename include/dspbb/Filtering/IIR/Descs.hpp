#pragma once

#include <type_traits>
#include <utility>


namespace dspbb {


namespace impl {
	struct IirMethod {};
	struct IirMethodButterworth : IirMethod {};


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
		ParamType low = ParamType(0.25);
		ParamType high = ParamType(0.75);

		template <class NewParamType>
		[[nodiscard]] auto Band(NewParamType lowNew, NewParamType highNew) {
			return Desc<IirMethodButterworth, NewParamType>{ { std::move(lowNew), std::move(highNew) } };
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

} // namespace impl


constexpr impl::IirMethodButterworth BUTTERWORTH;

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