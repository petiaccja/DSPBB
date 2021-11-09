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

template <class T, size_t Order, class TransformFunc>
FactoredPolynomial<T> TransformRoots(const FactoredPolynomial<T>& poly, TransformFunc func, size_t numRoots = 0, std::array<std::complex<T>, Order> padRoots = {}) {
	const auto [realToReals, realToPairs] = impl::CountTransformedRoots(poly.RealRoots().begin(), poly.RealRoots().end(), func);
	const auto [padReals, padPairs] = impl::CountRoots(padRoots.begin(), padRoots.end());
	if (numRoots >0  && numRoots < poly.NumRoots()) {
		throw std::invalid_argument("Number of transformed roots must be >= than number of input roots, or 0 to automatically deduce number.");
	}
	if (numRoots == 0) {
		numRoots = poly.NumRoots();
	}
	const size_t numPadSets = numRoots - poly.NumRoots();

	const size_t numReals = realToReals + numPadSets * padReals;
	const size_t numComplexPairs = realToPairs + Order * poly.NumComplexPairs() + numPadSets * padPairs;

	FactoredPolynomial<T> r;
	r.Resize(numReals, numComplexPairs);
	auto realRootIt = r.RealRoots().begin();
	auto complexRootIt = r.ComplexPairs().begin();

	for (auto& root : poly.RealRoots()) {
		const auto transformedRoots = func(root);
		assert(std::distance(std::begin(transformedRoots), std::end(transformedRoots)) == Order);
		for (const auto& troot : transformedRoots) {
			if (imag(troot) > T(0)) {
				*(complexRootIt++) = troot;
			}
			else if (imag(troot) == T(0)) {
				*(realRootIt++) = real(troot);
			}
		}
	}
	for (auto& root : poly.ComplexPairs()) {
		const auto transformedRoots = func(root);
		assert(std::distance(std::begin(transformedRoots), std::end(transformedRoots)) == Order);
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
		   * std::transform_reduce(poly.ComplexPairs().begin(), poly.ComplexPairs().end(), T(1), std::multiplies<>{}, pairGain);
}

} // namespace dspbb