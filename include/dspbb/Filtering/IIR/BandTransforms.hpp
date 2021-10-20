#pragma once

#include "../../LTISystems/System.hpp"

namespace dspbb {

template <class T>
PoleZeroSystem<T, eSystemDiscretization::CONTINUOUS> ScaleFrequency(const PoleZeroSystem<T, eSystemDiscretization::CONTINUOUS>& system, T scale) {
	auto zerosScaled = system.zeros;
	auto polesScaled = system.poles;
	zerosScaled.RealRoots() *= scale;
	zerosScaled.ComplexRoots() *= scale;
	polesScaled.RealRoots() *= scale;
	polesScaled.ComplexRoots() *= scale;

	const T gain = std::pow(scale, T(polesScaled.NumRoots() - zerosScaled.NumRoots()));

	return { system.gain * gain, zerosScaled, polesScaled };
}

} // namespace dspbb