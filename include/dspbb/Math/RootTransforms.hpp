#pragma once

#include "Polynomials.hpp"

#include <complex>

namespace dspbb {

namespace impl {

	template <class Iter>
	std::pair<size_t, size_t> CountRoots(Iter first, Iter last) {
		using T = remove_complex_t<std::decay_t<decltype(*first)>>;
		size_t realToReals = 0;
		size_t realToPairs = 0;
		assert(impl::AreRootsConjugatePairs(first, last)); // Mistake within DSPBB's transform functions, no exception.
		for (auto it = first; it != last; ++it) {
			const auto& root = *it;
			realToReals += size_t(imag(root) == T(0));
			realToPairs += size_t(imag(root) > T(0));
		}
		return { realToReals, realToPairs };
	}

	template <class Iter, class TransformFunc>
	std::pair<size_t, size_t> CountTransformedRoots(Iter first, Iter last, TransformFunc func) {
		size_t realToReals = 0;
		size_t realToPairs = 0;
		for (auto it = first; it != last; ++it) {
			const auto& root = *it;
			const auto transformedRoots = func(root);
			const auto [partialRealToReals, partialRealToPairs] = CountRoots(transformedRoots.begin(), transformedRoots.end());
			realToReals += partialRealToReals;
			realToPairs += partialRealToPairs;
		}
		return { realToReals, realToPairs };
	}

} // namespace impl

template <class T, size_t Order, class TransformFunc>
FactoredPolynomial<T> TransformRoots(const FactoredPolynomial<T>& poly, TransformFunc func, size_t numRoots = 0, std::array<std::complex<T>, Order> padRoots = {}) {
	const auto [realToReals, realToPairs] = impl::CountTransformedRoots(poly.RealRoots().begin(), poly.RealRoots().end(), func);
	const auto [padReals, padPairs] = impl::CountRoots(padRoots.begin(), padRoots.end());
	const size_t numPadSets = numRoots - poly.NumRoots();

	const size_t numReals = realToReals + numPadSets * padReals;
	const size_t numComplexPairs = realToPairs + Order * poly.NumComplexPairs() + numPadSets * padPairs;

	FactoredPolynomial<T> r;
	r.Resize(numReals, numComplexPairs);
	auto realRootIt = r.RealRoots().begin();
	auto complexRootIt = r.ComplexRoots().begin();

	for (auto& root : poly.RealRoots()) {
		const auto transformedRoots = func(root);
		for (const auto& troot : transformedRoots) {
			if (imag(troot) > T(0)) {
				*(complexRootIt++) = troot;
			}
			else if (imag(troot) == T(0)) {
				*(realRootIt++) = real(troot);
			}
		}
	}
	for (auto& root : poly.ComplexRoots()) {
		const auto transformedRoots = func(root);
		complexRootIt = std::copy(transformedRoots.begin(), transformedRoots.end(), complexRootIt);
	}
	for (size_t i = 0; i < numPadSets; ++i) {
		for (const auto& proot : padRoots) {
			if (imag(proot) > T(0)) {
				*(complexRootIt++) = proot;
			}
			else if (imag(proot) == T(0)) {
				*(realRootIt++) = real(proot);
			}
		}
	}

	return r;
}

template <class T, class RealGain, class PairGain>
T TransformGain(const FactoredPolynomial<T>& poly, RealGain realGain, PairGain pairGain) {
	return std::transform_reduce(poly.RealRoots().begin(), poly.RealRoots().end(), T(1), std::multiplies<>{}, realGain)
		   * std::transform_reduce(poly.ComplexRoots().begin(), poly.ComplexRoots().end(), T(1), std::multiplies<>{}, pairGain);
}

} // namespace dspbb