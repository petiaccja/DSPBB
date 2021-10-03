#pragma once

#include "../../LTISystems/PoleZeroSystem.hpp"

namespace dspbb {

template <class T>
PoleZeroSystem<T, eSystemDiscretization::CONTINUOUS> ScaleFrequency(const PoleZeroSystem<T, eSystemDiscretization::CONTINUOUS>& system, T scale) {
	auto polesScaled = system.Poles();
	auto zerosScaled = system.Zeros();
	T gain = T(1.0f);
	for (auto& p : polesScaled) {
		p *= scale;
		gain *= scale;
	}
	for (auto& z : zerosScaled) {
		z *= scale;
		gain /= scale;
	}
	return { system.Gain() * gain, polesScaled, zerosScaled };
}

} // namespace dspbb