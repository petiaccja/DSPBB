#pragma once

#include "PoleZeroSystem.hpp"

#include <optional>
#include <stdexcept>

namespace dspbb {


template <class T>
ContinuousPoleZeroSystem<T> BilinearTransform(const DiscretePoleZeroSystem<T>& discrete, T sampleRate) {
	throw std::logic_error{ "not implemented" };
}

template <class T>
DiscretePoleZeroSystem<T> BilinearTransform(const ContinuousPoleZeroSystem<T>& continuous, T sampleRate, std::optional<T> prewarp = {}) {
	const T k = prewarp ? *prewarp / std::tan(*prewarp * sampleRate / T(2)) : T(2) * sampleRate;

	const T gain = continuous(k);
	std::vector<std::complex<T>> poles;
	std::vector<std::complex<T>> zeros;
	poles.reserve(continuous.Poles().size());
	zeros.reserve(poles.capacity());

	// TODO: use ranges transform.
	std::transform(continuous.Poles().begin(), continuous.Poles().end(), std::back_inserter(poles), [&k](const auto& pole) {
		return (k + pole) / (k - pole);
	});
	std::transform(continuous.Zeros().begin(), continuous.Zeros().end(), std::back_inserter(zeros), [&k](const auto& zero) {
		return (k + zero) / (k - zero);
	});
	zeros.resize(poles.size(), -T(1));

	return { gain, poles, zeros };
}


} // namespace dspbb