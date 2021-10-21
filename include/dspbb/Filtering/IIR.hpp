#pragma once

#include "../LTISystems/DiscretizationTransforms.hpp"
#include "../LTISystems/System.hpp"
#include "../LTISystems/RepresentationTransforms.hpp"
#include "IIR/BandTransforms.hpp"
#include "IIR/Butterworth.hpp"
#include "IIR/Descs.hpp"
#include "IIR/Filter.hpp"
#include "IIR/Realizations.hpp"


namespace dspbb {


//------------------------------------------------------------------------------
// Butterworth method
//------------------------------------------------------------------------------

// Lowpass
template <class T, class ParamType>
auto IirFilter(size_t order, const impl::LowpassDesc<impl::IirMethodButterworth, ParamType>& desc) {
	constexpr T sampleRate = T(1);
	constexpr T angularLimit = sampleRate * pi_v<T>;
	const auto system = Butterworth<ParamType>(order);
	const auto scaled = ScaleFrequency(system, pi_v<T> * desc.cutoff);
	auto discrete = BilinearTransform(scaled, sampleRate, { T(desc.cutoff) * angularLimit });
	return discrete;
}

template <class T, class ParamType>
void IirFilter(DiscretePoleZeroSystem<T>& out, const impl::LowpassDesc<impl::IirMethodButterworth, ParamType>& desc) {
	assert(out.zeros.NumRoots() == out.poles.NumRoots());
	const size_t order = out.poles.NumRoots();
}

} // namespace dspbb