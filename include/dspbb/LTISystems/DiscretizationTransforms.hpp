#pragma once

#include "Systems.hpp"
#include "../Math/RootTransforms.hpp"

#include <optional>
#include <stdexcept>

namespace dspbb {


template <class T>
ContinuousZeroPoleGain<T> BilinearTransform(const DiscreteZeroPoleGain<T>& discrete, T sampleRate) {
	throw std::logic_error{ "not implemented" };
}

template <class T>
DiscreteZeroPoleGain<T> BilinearTransform(const ContinuousZeroPoleGain<T>& continuous, T sampleRate, std::optional<T> prewarp = {}) {
	const T k = prewarp ? *prewarp / std::tan(*prewarp / sampleRate / T(2)) : T(2) * sampleRate;
	
	// Opposite roots
	const std::complex<T> z = -T(1);
	const std::array z12 = { z };

	// Roots
	const auto transform = [k](const auto& p) {
		std::complex<T> z = (k + p) / (k - p);
		return std::array{ z };
	};

	// Do transform
	const size_t numRoots = std::max(continuous.zeros.NumRoots(), continuous.poles.NumRoots());
	FactoredPolynomial<T> newZeros = TransformRoots(continuous.zeros, transform, numRoots, z12);
	FactoredPolynomial<T> newPoles = TransformRoots(continuous.poles, transform, numRoots, z12);
	const T newGain = continuous(k);

	return { newGain, std::move(newZeros), std::move(newPoles) };
}


} // namespace dspbb