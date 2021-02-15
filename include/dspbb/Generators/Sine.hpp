#pragma once

#include <cstdint>
#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Utility/Numbers.hpp>
#include <dspbb/Utility/TypeTraits.hpp>


namespace dspbb {

template <class T, eSignalDomain Domain>
Signal<T, Domain> SineWave(size_t length, uint64_t sampleRate, remove_complex_t<T> frequency, remove_complex_t<T> phase = 0) {
	Signal<T, Domain> signal(length);
	for (size_t i = 0; i < length; ++i) {
		double time = double(i) / double(sampleRate);
		double totalPhase = 2.0 * pi_v<double> * time * double(frequency) + double(phase);
		signal[i] = T(remove_complex_t<T>(std::sin(totalPhase)));
	}
	return signal;
}



} // namespace dspbb
