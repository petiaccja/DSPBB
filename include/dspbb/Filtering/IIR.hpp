#pragma once

#include "../LTISystems/DiscretizationTransforms.hpp"
#include "../LTISystems/Systems.hpp"
#include "FilterUtility.hpp"
#include "IIR/BandTransforms.hpp"
#include "IIR/Butterworth.hpp"
#include "IIR/Chebyshev.hpp"
#include "IIR/Descs.hpp"
#include "IIR/Elliptic.hpp"
#include "IIR/Filter.hpp"
#include "IIR/Realizations.hpp"


namespace dspbb {


//------------------------------------------------------------------------------
// Prototype filters
//------------------------------------------------------------------------------

namespace impl {
	template <class T>
	DiscreteZeroPoleGain<T> DiscretizePrototype(const ContinuousZeroPoleGain<T>& sys) {
		constexpr T sampleRate = T(2) / pi_v<T>;
		auto halfband = BilinearTransform(sys, sampleRate, { T(1) });
		return halfband;
	}

	template <class T>
	auto PrototypeButterworth(size_t order) {
		return DiscretizePrototype(Butterworth<T>(order));
	}

	template <class T>
	auto PrototypeChebyshev1(size_t order, T passbandRipple) {
		return DiscretizePrototype(Chebyshev1<T>(order, passbandRipple));
	}

	template <class T>
	auto PrototypeChebyshev2(size_t order, T stopbandRipple) {
		return DiscretizePrototype(Chebyshev2<T>(order, stopbandRipple));
	}

	template <class T>
	auto PrototypeElliptic(size_t order, T passbandRipple, T stopbandRipple) {
		return DiscretizePrototype(Elliptic<T>(order, passbandRipple, stopbandRipple));
	}
} // namespace impl


//------------------------------------------------------------------------------
// Butterworth method
//------------------------------------------------------------------------------

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::butterworth::LowpassDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeButterworth<T>(order);
	auto filter = Halfband2Lowpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::butterworth::HighpassDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeButterworth<T>(order);
	auto filter = Halfband2Highpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::butterworth::BandpassDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeButterworth<T>(order / 2);
	auto filter = Halfband2Bandpass(halfband, desc.lower, desc.upper);
	return filter;
}

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::butterworth::BandstopDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeButterworth<T>(order / 2);
	auto filter = Halfband2Bandstop(halfband, desc.lower, desc.upper);
	return filter;
}

//------------------------------------------------------------------------------
// Chebyshev 1 method
//------------------------------------------------------------------------------

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::chebyshev1::LowpassDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeChebyshev1<T>(order, desc.passbandRipple);
	auto filter = Halfband2Lowpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::chebyshev1::HighpassDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeChebyshev1<T>(order, desc.passbandRipple);
	auto filter = Halfband2Highpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::chebyshev1::BandpassDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeChebyshev1<T>(order / 2, desc.passbandRipple);
	auto filter = Halfband2Bandpass(halfband, desc.lower, desc.upper);
	return filter;
}

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::chebyshev1::BandstopDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeChebyshev1<T>(order / 2, desc.passbandRipple);
	auto filter = Halfband2Bandstop(halfband, desc.lower, desc.upper);
	return filter;
}

//------------------------------------------------------------------------------
// Chebyshev 2 method
//------------------------------------------------------------------------------

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::chebyshev2::LowpassDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeChebyshev2<T>(order, desc.stopbandRipple);
	auto filter = Halfband2Lowpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::chebyshev2::HighpassDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeChebyshev2<T>(order, desc.stopbandRipple);
	auto filter = Halfband2Highpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::chebyshev2::BandpassDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeChebyshev2<T>(order / 2, desc.stopbandRipple);
	auto filter = Halfband2Bandpass(halfband, desc.lower, desc.upper);
	return filter;
}

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::chebyshev2::BandstopDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeChebyshev2<T>(order / 2, desc.stopbandRipple);
	auto filter = Halfband2Bandstop(halfband, desc.lower, desc.upper);
	return filter;
}

//------------------------------------------------------------------------------
// Elliptic method
//------------------------------------------------------------------------------

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::elliptic::LowpassDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeElliptic<T>(order, desc.passbandRipple, desc.stopbandRipple);
	auto filter = Halfband2Lowpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::elliptic::HighpassDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeElliptic<T>(order, desc.passbandRipple, desc.stopbandRipple);
	auto filter = Halfband2Highpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::elliptic::BandpassDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeElliptic<T>(order / 2, desc.passbandRipple, desc.stopbandRipple);
	auto filter = Halfband2Bandpass(halfband, desc.lower, desc.upper);
	return filter;
}

template <class T, class ParamType>
auto DesignFilter(size_t order, const impl::elliptic::BandstopDesc<ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeElliptic<T>(order / 2, desc.passbandRipple, desc.stopbandRipple);
	auto filter = Halfband2Bandstop(halfband, desc.lower, desc.upper);
	return filter;
}

//------------------------------------------------------------------------------
// In-place wrapper
//------------------------------------------------------------------------------

template <class T, class ResponseDesc>
void DesignFilter(DiscreteZeroPoleGain<T>& out, const ResponseDesc& desc) {
	// TODO:
	// This is kinda useless. I should maybe add view-like analogs to LTI systems.
	assert(out.zeros.num_roots() == out.poles.num_roots());
	const size_t order = out.poles.num_roots();
	const auto filter = DesignFilter(order, desc);
	out.regroup(filter.num_real_roots());
	std::copy(filter.real_roots().begin(), filter.real_roots.end(), out.real_roots.begin());
	std::copy(filter.complex_pairs().begin(), filter.complex_pairs.end(), out.complex_pairs.begin());
}

} // namespace dspbb