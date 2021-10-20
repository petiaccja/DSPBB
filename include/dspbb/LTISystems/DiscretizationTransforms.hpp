#pragma once

#include "System.hpp"

#include <optional>
#include <stdexcept>

namespace dspbb {


template <class T>
ContinuousPoleZeroSystem<T> BilinearTransform(const DiscretePoleZeroSystem<T>& discrete, T sampleRate) {
	throw std::logic_error{ "not implemented" };
}

template <class T>
DiscretePoleZeroSystem<T> BilinearTransform(const ContinuousPoleZeroSystem<T>& continuous, T sampleRate, std::optional<T> prewarp = {}) {
	const T k = prewarp ? *prewarp / std::tan(*prewarp / sampleRate / T(2)) : T(2) * sampleRate;

	const T gain = continuous(k);
	FactoredPolynomial<T> poles;
	FactoredPolynomial<T> zeros;

	const size_t numDiscreteRoots = std::max(continuous.zeros.NumRoots(), continuous.poles.NumRoots());
	const size_t numZeroPairs = continuous.zeros.NumComplexPairs();
	const size_t numPolePairs = continuous.poles.NumComplexPairs();
	zeros.Resize(numDiscreteRoots - numZeroPairs * 2, numZeroPairs, -T(1));
	poles.Resize(numDiscreteRoots - numPolePairs * 2, numPolePairs, -T(1));

	// TODO: use ranges transform.
	std::transform(continuous.zeros.RealRoots().begin(), continuous.zeros.RealRoots().end(), zeros.RealRoots().begin(), [&k](const auto& zero) {
		return (k + zero) / (k - zero);
	});
	std::transform(continuous.zeros.ComplexRoots().begin(), continuous.zeros.ComplexRoots().end(), zeros.ComplexRoots().begin(), [&k](const auto& zero) {
		return (k + zero) / (k - zero);
	});
	std::transform(continuous.poles.RealRoots().begin(), continuous.poles.RealRoots().end(), poles.RealRoots().begin(), [&k](const auto& pole) {
		return (k + pole) / (k - pole);
	});
	std::transform(continuous.poles.ComplexRoots().begin(), continuous.poles.ComplexRoots().end(), poles.ComplexRoots().begin(), [&k](const auto& pole) {
		return (k + pole) / (k - pole);
	});

	return { gain, std::move(zeros), std::move(poles) };
}


} // namespace dspbb