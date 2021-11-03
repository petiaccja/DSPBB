#pragma once

#include "../../LTISystems/System.hpp"
#include "../../Utility/Numbers.hpp"
#include "../../Math/RootTransforms.hpp"

#include <complex>

namespace dspbb {

namespace impl {
	
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
		FactoredPolynomial<T> newZeros = TransformRoots(system.zeros, transform, numRoots, z12);
		FactoredPolynomial<T> newPoles = TransformRoots(system.poles, transform, numRoots, z12);
		const T newGain = system.gain * TransformGain(system.zeros, rgain, cgain) / TransformGain(system.poles, rgain, cgain);

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
		FactoredPolynomial<T> newZeros = TransformRoots(system.zeros, transform, numRoots, z12);
		FactoredPolynomial<T> newPoles = TransformRoots(system.poles, transform, numRoots, z12);
		const T newGain = system.gain * TransformGain(system.zeros, rgain, cgain) / TransformGain(system.poles, rgain, cgain);

		return { newGain, std::move(newZeros), std::move(newPoles) };
	}

}

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