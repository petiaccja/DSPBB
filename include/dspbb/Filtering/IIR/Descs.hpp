#pragma once

#include <type_traits>


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

	template <class ParamType>
	struct LowpassDesc<IirMethodButterworth, ParamType> {
		ParamType cutoff;

		[[nodiscard]] auto Cutoff(ParamType cutoffNew) {
			return LowpassDesc<IirMethodButterworth, ParamType>{ cutoffNew };
		}
	};

	template <>
	struct LowpassDesc<IirMethodButterworth> {
		template <class ParamType>
		[[nodiscard]] auto Cutoff(ParamType cutoffNew) {
			return LowpassDesc<IirMethodButterworth, ParamType>{ cutoffNew };
		}
	};

} // namespace impl


constexpr impl::IirMethodButterworth BUTTERWORTH;

template <class Method, std::enable_if_t<std::is_base_of_v<impl::IirMethod, Method>, int> = 0>
auto Lowpass(Method) {
	return impl::LowpassDesc<Method>{};
}


} // namespace dspbb