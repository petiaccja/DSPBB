#pragma once

#include "../../LTISystems/System.hpp"
#include "../../Utility/Numbers.hpp"

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

	template <class T, size_t Order, class TransformFunc>
	FactoredPolynomial<T> TransformRoots(const FactoredPolynomial<T>& poly, TransformFunc func, size_t numRoots = 0, std::array<std::complex<T>, Order> padRoots = {}) {
		const auto [realToReals, realToPairs] = CountTransformedRoots(poly.RealRoots().begin(), poly.RealRoots().end(), func);
		const auto [padReals, padPairs] = CountRoots(padRoots.begin(), padRoots.end());
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


	template <class T>
	std::complex<T> MapZDomain(std::complex<T> z, T s, T a0, T a1) {
		return s * (a1 * z + a0) / (a1 + a0 * z);
	}

	template <class T>
	std::complex<T> MapZDomain(std::complex<T> z, T s, T a0, T a1, T a2) {
		return s * (a2 * z * z + a1 * z + a0) / (a2 + a1 * z + a0 * z * z);
	}

	template <class T>
	DiscretePoleZeroSystem<T> MapZDomain(const DiscretePoleZeroSystem<T>& system, T s, T a0, T a1) {
		// Opposite roots
		const std::complex<T> z = -(a1 / a0);
		const std::array z12 = { z };

		// Roots
		const auto transform = [s, a1, a0](const auto& p) {
			std::complex<T> z = (-(a1 * p) + a0 * s) / (a0 * p - a1 * s);
			return std::array{ z };
		};

		// Gain
		const auto cgain = [s, a1, a0](const std::complex<T>& p) {
			return real(((a1 * s) / a0 - p) * ((a1 * s) / a0 - conj(p)));
		};
		const auto rgain = [s, a1, a0](const T& p) {
			return (a1 * s) / a0 - p;
		};

		// Do transform
		const size_t numRoots = std::max(system.zeros.NumRoots(), system.poles.NumRoots());
		FactoredPolynomial<T> newZeros = impl::TransformRoots(system.zeros, transform, numRoots, z12);
		FactoredPolynomial<T> newPoles = impl::TransformRoots(system.poles, transform, numRoots, z12);
		const T newGain = system.gain * impl::TransformGain(system.zeros, rgain, cgain) / impl::TransformGain(system.poles, rgain, cgain);

		return { newGain, std::move(newZeros), std::move(newPoles) };
	}

	template <class T>
	DiscretePoleZeroSystem<T> MapZDomain(const DiscretePoleZeroSystem<T>& system, T s, T a0, T a1, T a2) {
		// Opposite roots
		const std::complex<T> det = a1 * a1 - T(4) * a0 * a2;
		const std::complex<T> z1 = (-a1 - std::sqrt(det)) / (T(2) * a0);
		const std::complex<T> z2 = (-a1 + std::sqrt(det)) / (T(2) * a0);
		const std::array z12 = { z1, z2 };

		// Roots
		const auto transform = [s, a2, a1, a0](const auto& p) {
			const std::complex<T> det = (a1 * p - a1 * s) * (a1 * p - a1 * s) + T(4) * (a2 * p - a0 * s) * (-(a0 * p) + a2 * s);
			std::complex<T> z1 = -((a1 * p - a1 * s + std::sqrt(det)) / (2 * a0 * p - 2 * a2 * s));
			std::complex<T> z2 = -((a1 * p - a1 * s - std::sqrt(det)) / (2 * a0 * p - 2 * a2 * s));
			return std::array{ z1, z2 };
		};

		// Gain
		const auto cgain = [s, a2, a0](const std::complex<T>& p) {
			return real((a2 * s / a0 - p) * (a2 * s / a0 - conj(p)));
		};
		const auto rgain = [s, a2, a0](const T& p) {
			return a2 * s / a0 - p;
		};

		// Do transform
		const size_t numRoots = std::max(system.zeros.NumRoots(), system.poles.NumRoots());
		FactoredPolynomial<T> newZeros = impl::TransformRoots(system.zeros, transform, numRoots, z12);
		FactoredPolynomial<T> newPoles = impl::TransformRoots(system.poles, transform, numRoots, z12);
		const T newGain = system.gain * impl::TransformGain(system.zeros, rgain, cgain) / impl::TransformGain(system.poles, rgain, cgain);

		return { newGain, std::move(newZeros), std::move(newPoles) };
	}

} // namespace impl


template <class T>
DiscretePoleZeroSystem<T> Halfband2Lowpass(const DiscretePoleZeroSystem<T>& system, T to) {
	const T w = to * pi_v<T>;

	const T s = T(1);
	const T a1 = 1;
	const T a0 = -(std::cos(w) / (1 + std::sin(w)));

	return impl::MapZDomain(system, s, a0, a1);
}

template <class T>
DiscretePoleZeroSystem<T> Halfband2Highpass(const DiscretePoleZeroSystem<T>& system, T to) {
	const T w = to * pi_v<T>;

	const T s = -T(1);
	const T a1 = std::cos(w) / (-1 + std::sin(w));
	const T a0 = 1;

	return impl::MapZDomain(system, s, a0, a1);
}

template <class T>
DiscretePoleZeroSystem<T> Halfband2Bandpass(const DiscretePoleZeroSystem<T>& system, T to1, T to2) {
	const T w1 = to1 * pi_v<T>;
	const T w2 = to2 * pi_v<T>;

	const T s = -T(1);
	const T a2 = -1 + 2 / (1 + std::tan((w1 - w2) / 2));
	const T a1 = -((std::cos(w1) + std::cos(w2) + std::sin(w1) - std::sin(w2)) / (1 + std::sin(w1 - w2)));
	const T a0 = 1;

	return impl::MapZDomain(system, s, a0, a1, a2);
}

template <class T>
DiscretePoleZeroSystem<T> Halfband2Bandstop(const DiscretePoleZeroSystem<T>& system, T to1, T to2) {
	const T w1 = to1 * pi_v<T>;
	const T w2 = to2 * pi_v<T>;

	const T s = T(1);
	const T a2 = 1;
	const T a1 = (std::cos(w1) + std::cos(w2) - std::sin(w1) + std::sin(w2)) / (-1 + std::sin(w1 - w2));
	const T a0 = -1 - 2 / (-1 + std::tan((w1 - w2) / 2));

	return impl::MapZDomain(system, s, a0, a1, a2);
}


} // namespace dspbb