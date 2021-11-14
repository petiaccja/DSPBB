#pragma once

#include "../LTISystems/DiscretizationTransforms.hpp"
#include "../LTISystems/Systems.hpp"
#include "FilterUtility.hpp"
#include "IIR/BandTransforms.hpp"
#include "IIR/Butterworth.hpp"
#include "IIR/Descs.hpp"
#include "IIR/Filter.hpp"
#include "IIR/Realizations.hpp"


namespace dspbb {


//------------------------------------------------------------------------------
// Butterworth method
//------------------------------------------------------------------------------

namespace impl {
	template <class T>
	auto PrototypeButterworth(size_t order) {
		constexpr T sampleRate = T(2) / pi_v<T>;
		const auto analog = Butterworth<T>(order);
		auto halfband = BilinearTransform(analog, sampleRate, { T(1) });
		return halfband;
	}
} // namespace impl

// Lowpass
template <class T, class ParamType>
auto IirFilter(size_t order, const impl::LowpassDesc<impl::IirMethodButterworth, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeButterworth<T>(order);
	auto filter = Halfband2Lowpass(halfband, desc.cutoff);
	return filter;
}

// Highpass
template <class T, class ParamType>
auto IirFilter(size_t order, const impl::HighpassDesc<impl::IirMethodButterworth, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	const auto halfband = impl::PrototypeButterworth<T>(order);
	auto filter = Halfband2Highpass(halfband, desc.cutoff);
	return filter;
}

// Bandpass
template <class T, class ParamType>
auto IirFilter(size_t order, const impl::BandpassDesc<impl::IirMethodButterworth, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	if (order % 2 != 0) {
		throw std::invalid_argument("IIR bandpass filter must have an even order.");
	}
	const auto halfband = impl::PrototypeButterworth<T>(order);
	auto filter = Halfband2Bandpass(halfband, desc.low, desc.high);
	return filter;
}

// Highpass
template <class T, class ParamType>
auto IirFilter(size_t order, const impl::BandstopDesc<impl::IirMethodButterworth, ParamType>& desc) -> DiscreteZeroPoleGain<T> {
	if (order % 2 != 0) {
		throw std::invalid_argument("IIR bandstop filter must have an even order.");
	}
	const auto halfband = impl::PrototypeButterworth<T>(order);
	auto filter = Halfband2Bandstop(halfband, desc.low, desc.high);
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