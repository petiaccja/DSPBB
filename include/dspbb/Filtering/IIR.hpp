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
auto IirFilter(size_t order, const impl::LowpassDesc<impl::IirMethodButterworth, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeButterworth<T>(order);
	auto filter = Halfband2Lowpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::HighpassDesc<impl::IirMethodButterworth, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeButterworth<T>(order);
	auto filter = Halfband2Highpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::BandpassDesc<impl::IirMethodButterworth, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeButterworth<T>(order / 2);
	auto filter = Halfband2Bandpass(halfband, desc.lower, desc.upper);
	return filter;
}

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::BandstopDesc<impl::IirMethodButterworth, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeButterworth<T>(order / 2);
	auto filter = Halfband2Bandstop(halfband, desc.lower, desc.upper);
	return filter;
}

//------------------------------------------------------------------------------
// Chebyshev 1 method
//------------------------------------------------------------------------------

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::LowpassDesc<impl::IirMethodChebyshev1, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeChebyshev1<T>(order, desc.passbandRipple);
	auto filter = Halfband2Lowpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::HighpassDesc<impl::IirMethodChebyshev1, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeChebyshev1<T>(order, desc.passbandRipple);
	auto filter = Halfband2Highpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::BandpassDesc<impl::IirMethodChebyshev1, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeChebyshev1<T>(order / 2, desc.passbandRipple);
	auto filter = Halfband2Bandpass(halfband, desc.lower, desc.upper);
	return filter;
}

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::BandstopDesc<impl::IirMethodChebyshev1, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeChebyshev1<T>(order / 2, desc.passbandRipple);
	auto filter = Halfband2Bandstop(halfband, desc.lower, desc.upper);
	return filter;
}

//------------------------------------------------------------------------------
// Chebyshev 2 method
//------------------------------------------------------------------------------

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::LowpassDesc<impl::IirMethodChebyshev2, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeChebyshev2<T>(order, desc.stopbandRipple);
	auto filter = Halfband2Lowpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::HighpassDesc<impl::IirMethodChebyshev2, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeChebyshev2<T>(order, desc.stopbandRipple);
	auto filter = Halfband2Highpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::BandpassDesc<impl::IirMethodChebyshev2, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeChebyshev2<T>(order / 2, desc.stopbandRipple);
	auto filter = Halfband2Bandpass(halfband, desc.lower, desc.upper);
	return filter;
}

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::BandstopDesc<impl::IirMethodChebyshev2, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeChebyshev2<T>(order / 2, desc.stopbandRipple);
	auto filter = Halfband2Bandstop(halfband, desc.lower, desc.upper);
	return filter;
}

//------------------------------------------------------------------------------
// Elliptic method
//------------------------------------------------------------------------------

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::LowpassDesc<impl::IirMethodElliptic, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeElliptic<T>(order, desc.passbandRipple, desc.stopbandRipple);
	auto filter = Halfband2Lowpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::HighpassDesc<impl::IirMethodElliptic, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeElliptic<T>(order, desc.passbandRipple, desc.stopbandRipple);
	auto filter = Halfband2Highpass(halfband, desc.cutoff);
	return filter;
}

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::BandpassDesc<impl::IirMethodElliptic, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeElliptic<T>(order / 2, desc.passbandRipple, desc.stopbandRipple);
	auto filter = Halfband2Bandpass(halfband, desc.lower, desc.upper);
	return filter;
}

template <class T, class ParamType>
auto IirFilter(size_t order, const impl::BandstopDesc<impl::IirMethodElliptic, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	assert(order % 2 == 0);
	const auto halfband = impl::PrototypeElliptic<T>(order / 2, desc.passbandRipple, desc.stopbandRipple);
	auto filter = Halfband2Bandstop(halfband, desc.lower, desc.upper);
	return filter;
}

//------------------------------------------------------------------------------
// In-place wrapper
//------------------------------------------------------------------------------

template <class T, class ResponseDesc>
void IirFilter(DiscreteZeroPoleGain<T>& out, const ResponseDesc& desc) {
	// TODO:
	// This is kinda useless. I should maybe add view-like analogs to LTI systems.
	assert(out.zeros.NumRoots() == out.poles.NumRoots());
	const size_t order = out.poles.NumRoots();
	const auto filter = IirFilter(order, desc);
	out.Regroup(filter.NumRealRoots());
	std::copy(filter.RealRoots().begin(), filter.RealRoots.end(), out.RealRoots.begin());
	std::copy(filter.ComplexPairs().begin(), filter.ComplexPairs.end(), out.ComplexPairs.begin());
}

} // namespace dspbb