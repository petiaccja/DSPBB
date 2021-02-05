#pragma once

#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Utility/Numbers.hpp>

#include <cstdint>


namespace enl {

template <class T, eSignalDomain Domain>
Signal<T, Domain> SineWave(size_t length, uint64_t sampleRate, double frequency, double phase = 0) {
	Signal<T, Domain> signal(length);
	for (size_t i=0; i<length; ++i) {
		double time = double(i) / double(sampleRate);
		double totalPhase = 2 * pi_v<double> * time * frequency + phase;
		signal[i] = T(std::sin(totalPhase));
	}
	return signal;
}



}
