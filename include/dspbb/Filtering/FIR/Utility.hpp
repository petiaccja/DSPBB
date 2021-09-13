#pragma once


namespace dspbb {

template <class T, class U>
T NormalizedFrequency(T frequency, U sampleRate) {
	return T(2) * frequency / T(sampleRate);
}

} // namespace dspbb