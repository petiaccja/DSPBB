#pragma once

#include "../../Generators/Spaces.hpp"
#include "../../LTISystems/Systems.hpp"
#include "../../Math/EllipticFunctions.hpp"
#include "../../Math/Solvers.hpp"
#include "../../Utility/Numbers.hpp"


namespace dspbb {

namespace impl {

	template <class T>
	struct EllipticParameters {
		T k;
		T kp;
		T K;
		T Kp;
		T k1;
		T k1p;
		T K1;
		T K1p;
		T epsilon;
	};

	template <class T>
	EllipticParameters<T> EllipticOrderRipples(size_t order, T passbandRipple, T stopbandRipple) {
		const T epsilon = std::sqrt((T(2) * passbandRipple - passbandRipple * passbandRipple)
									/ (T(1) - T(2) * passbandRipple + passbandRipple * passbandRipple));
		const T k1 = epsilon / std::sqrt(T(1) / (stopbandRipple * stopbandRipple) - T(1));
		const T k1p = std::sqrt(T(1) - k1 * k1);
		const T K1 = EllipticK(k1);
		const T K1p = EllipticK(k1p);

		const T Kratio = T(order) * K1 / K1p;

		// Bisecting may not be too efficient...
		// There is an algorithm in A. H. Gray and J. D. Markel, "A Computer Program for Designing Digital Elliptic Filters"
		// But I have no damn clue how that works, or even if it does what I want
		const T k = Bisect([&Kratio](auto k) { return EllipticK(k) / EllipticK(std::sqrt(T(1) - k * k)) - Kratio; }, T(0), T(1));
		const T kp = std::sqrt(T(1) - k * k);
		const T K = EllipticK(k);
		const T Kp = EllipticK(kp);

		return { k, kp, K, Kp, k1, k1p, K1, K1p, epsilon };
	}

} // namespace impl


template <class T>
ZeroPoleGain<T, eDiscretization::CONTINUOUS> Elliptic(size_t order, T passbandRipple, T stopbandRipple) {
	const auto [k, kp, K, Kp, k1, k1p, K1, K1p, epsilon] = impl::EllipticOrderRipples(order, passbandRipple, stopbandRipple);

	FactoredPolynomial<T> zeros;
	FactoredPolynomial<T> poles;
	zeros.resize(0, order / 2);
	poles.resize(order % 2, order / 2, -T(1));

	size_t index = 0;
	for (auto& root : zeros.complex_pairs()) {
		size_t i = 2 * index + 1 + order % 2;
		root = i_v<T> / (k * EllipticSN(T(i) * K / T(order), k));
		++index;
	}

	index = 0;
	const auto jv0 = EllipticArcSN(i_v<T> / epsilon, k1) * (K / T(order) / K1);
	for (auto& root : poles.real_roots()) {
		size_t i = 2 * index + 1 - order % 2;
		const auto croot = i_v<T> * EllipticSN(K * T(i) / T(order) + jv0, k);
		assert(std::abs(imag(croot)) < 0.001 * std::abs(real(croot)));
		root = real(croot);
		++index;
	}
	for (auto& root : poles.complex_pairs()) {
		size_t i = 2 * index + 1 - order % 2;
		root = i_v<T> * EllipticSN(K * T(i) / T(order) + jv0, k);
		++index;
	}

	const T gain = poles(T(0)) / zeros(T(0)) * (order % 2 ? T(1) : T(1) - passbandRipple);

	return { gain, std::move(zeros), std::move(poles) };
}


} // namespace dspbb