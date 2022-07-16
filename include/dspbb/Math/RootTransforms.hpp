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
		if (!impl::AreRootsConjugatePairs(first, last)) {
			throw std::invalid_argument("Multiple roots must be reals or complex conjugate pairs.");
		}
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

template <class T, size_t order, class TransformFunc>
FactoredPolynomial<T> TransformRoots(const FactoredPolynomial<T>& poly, TransformFunc func, size_t numRoots = 0, std::array<std::complex<T>, order> padRoots = {}) {
	const auto [realToReals, realToPairs] = impl::CountTransformedRoots(poly.real_roots().begin(), poly.real_roots().end(), func);
	const auto [padReals, padPairs] = impl::CountRoots(padRoots.begin(), padRoots.end());
	assert(numRoots == 0 || numRoots >= poly.num_roots());
	if (numRoots == 0) {
		numRoots = poly.num_roots();
	}
	const size_t numPadSets = numRoots - poly.num_roots();

	const size_t numReals = realToReals + numPadSets * padReals;
	const size_t numComplexPairs = realToPairs + order * poly.num_complex_pairs() + numPadSets * padPairs;

	FactoredPolynomial<T> r;
	r.resize(numReals, numComplexPairs);
	auto realRootIt = r.real_roots().begin();
	auto complexRootIt = r.complex_pairs().begin();

	for (auto& root : poly.real_roots()) {
		const auto transformedRoots = func(root);
		assert(std::distance(std::begin(transformedRoots), std::end(transformedRoots)) == order);
		for (const auto& troot : transformedRoots) {
			if (imag(troot) > T(0)) {
				*(complexRootIt++) = troot;
			}
			else if (imag(troot) == T(0)) {
				*(realRootIt++) = real(troot);
			}
		}
	}
	for (auto& root : poly.complex_pairs()) {
		const auto transformedRoots = func(root);
		assert(std::distance(std::begin(transformedRoots), std::end(transformedRoots)) == order);
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
	return std::transform_reduce(poly.real_roots().begin(), poly.real_roots().end(), T(1), std::multiplies<>{}, realGain)
		   * std::transform_reduce(poly.complex_pairs().begin(), poly.complex_pairs().end(), T(1), std::multiplies<>{}, pairGain);
}

} // namespace dspbb