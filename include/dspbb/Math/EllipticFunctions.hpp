#pragma once


#include "../Utility/Numbers.hpp"
#include "../Utility/TypeTraits.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <limits>


namespace dspbb {

//------------------------------------------------------------------------------
// Reference
//------------------------------------------------------------------------------

// Algorithms were gathered from NIST's DLMF:
// Elliptic integrals: https://dlmf.nist.gov/19
// Jacobi elliptic functions: https://dlmf.nist.gov/22



//------------------------------------------------------------------------------
// Carlson symmetric forms
//------------------------------------------------------------------------------

template <class T>
T CarlsonRF(T x, T y, T z) {
	using S = remove_complex_t<T>;
	constexpr S r = std::numeric_limits<S>::epsilon() / S(4);

	if (int(x == S(0)) + int(y == S(0)) + int(z == S(0)) >= 2) {
		return std::numeric_limits<S>::infinity();
	}

	T Amn = (x + y + z) / S(3);
	T An = Amn;
	S Qn = std::pow(S(3) * r, -S(1) / S(8)) * std::max({ std::abs(An - x), std::abs(An - y), std::abs(An - z) });
	T lambdan = T(0);

	while (Qn > std::abs(An)) {
		const auto lambdan2 = lambdan * lambdan;
		const auto l1 = std::sqrt(x * y + (x + y) * lambdan + lambdan2);
		const auto l2 = std::sqrt(x * z + (x + z) * lambdan + lambdan2);
		const auto l3 = std::sqrt(z * y + (z + y) * lambdan + lambdan2);
		lambdan += l1 + l2 + l3;
		x /= S(4);
		y /= S(4);
		z /= S(4);
		lambdan /= S(4);
		Qn /= S(4);
		Amn /= S(4);
		An = Amn + lambdan;
	}

	const auto X = (Amn - x) / An;
	const auto Y = (Amn - y) / An;
	const auto Z = (Amn - z) / An;
	const auto E2 = X * Y + X * Z + Y * Z;
	const auto E3 = X * Y * Z;
	const auto series = E2 * E2 * E3 / S(16)
						- S(5) * E2 * E2 * E2 / S(208)
						+ S(3) * E3 * E3 / S(104)
						- S(3) * E2 * E3 / S(44)
						+ E2 * E2 / S(24)
						+ E3 / S(14)
						- E2 / S(10)
						+ S(1);
	return series / std::sqrt(An);
}

template <class T1, class T2, class T3>
auto CarlsonRF(T1 x, T2 y, T3 z) -> std::common_type_t<T1, T2, T3> {
	using R = std::common_type_t<T1, T2, T3>;
	return CarlsonRF(R(x), R(y), R(z));
}


//------------------------------------------------------------------------------
// Elliptic integrals
//------------------------------------------------------------------------------

template <class T>
T EllipticK(T k) {
	return CarlsonRF(T(0), T(1) - k * k, T(1));
}


//------------------------------------------------------------------------------
// Jacobi amplitude function
//------------------------------------------------------------------------------

template <class T, std::enable_if_t<!is_complex_v<T>, int> = 0>
auto EllipticAM(T x, T k) {
	assert(T(0) <= k && k <= T(1));

	// Special cases
	if (k == T(1)) {
		return T(2) * std::atan(std::exp(x)) - pi_v<T> / 2;
	}
	if (k == T(0)) {
		return x;
	}

	// Iteration related constructs
	constexpr auto epsilon = std::numeric_limits<T>::epsilon() / T(4);
	constexpr auto iterationLimit = 32;
	thread_local std::array<T, iterationLimit> factors;

	// AGM forward iteration
	T an = 1;
	T bn = std::sqrt(T(1) - k * k);
	T cn = T(0.5) * (an - bn);
	int n = 0;

	while (std::abs(cn) > epsilon) {
		std::tie(an, bn, cn) = std::make_tuple(T(0.5) * (an + bn), std::sqrt(an * bn), T(0.5) * (an - bn));
		factors[n] = cn / an;
		++n;

		if (n >= iterationLimit) {
			throw std::domain_error("Arithmetic-geometric mean failed to converge within a reasonable number of iterations.");
		}
	}

	// Phi backwards iteration
	T phi = T(1 << n) * an * x;
	while (n-- > 0) {
		phi = T(0.5) * (phi + std::asin(factors[n] * std::sin(phi)));
	}

	return phi;
}

template <class T1, class T2, std::enable_if_t<!is_complex_v<T1> && !is_complex_v<T2>, int> = 0>
auto EllipticAM(T1 x, T2 k) -> std::common_type_t<T1, T2> {
	using R = std::common_type_t<T1, T2>;
	return EllipticAM(R(x), R(k));
}


//------------------------------------------------------------------------------
// Jacobi elliptic functions
//------------------------------------------------------------------------------


template <class T>
std::tuple<T, T, T, T, T> EllipticReduceRange(const T x, const T k) {
	T xPrime = x;
	T kPrime = k;
	T factorSn = T(1);
	T factorCn = T(1);
	T factorDn = T(1);

	if (kPrime < T(0)) {
		kPrime = -kPrime;
	}
	if (kPrime > T(1)) {
		xPrime *= kPrime;
		factorSn /= kPrime;
		kPrime = T(1) / kPrime;
	}

	if (xPrime > pi_v<T> / 2) {
		const auto K = EllipticK(k);
		int octant = 0;
		xPrime = std::remquo(x, T(2) * K, &octant);
		if (octant % 2 == 1) {
			factorSn *= T(-1);
			factorCn *= T(-1);
		}
	}

	return { xPrime, kPrime, factorSn, factorCn, factorDn };
}

template <class T, std::enable_if_t<!is_complex_v<T>, int> = 0>
std::tuple<T, T, T> EllipticSNCNDN(T x, T k) {
	const auto [xPrime, kPrime, factorSn, factorCn, factorDn] = EllipticReduceRange(x, k);

	const T am = EllipticAM(xPrime, kPrime);
	const T sn = factorSn * std::sin(am);
	const T cn = factorCn * std::cos(am);
	const T dn = factorDn * std::sqrt(T(1) - k * k * sn * sn);
	return { sn, cn, dn };
}

template <class T>
auto EllipticSNCNDN(std::complex<T> x, T k) -> std::tuple<std::complex<T>, std::complex<T>, std::complex<T>> {
	const auto [snr, cnr, dnr] = EllipticSNCNDN(real(x), k);
	const auto [sni, cni, dni] = EllipticSNCNDN(imag(x), std::sqrt(1 - k * k));

	const auto d = cni * cni + k * k * sni * sni * snr * snr;
	const auto sn = std::complex<T>(snr * dni, sni * cni * cnr * dnr) / d;
	const auto cn = std::complex<T>(cnr * cni, -snr * dnr * sni * dni) / d;
	const auto dn = std::complex<T>(dnr * dni * cni, -k * k * snr * cnr * sni) / d;
	return { sn, cn, dn };
}

template <class T1, class T2, class T = std::common_type_t<T1, T2>, std::enable_if_t<!is_complex_v<T1> && !is_complex_v<T2>, int> = 0>
auto EllipticSNCNDN(T1 x, T2 k) -> std::tuple<T, T, T> {
	return EllipticSNCNDN(T(x), T(k));
}

template <class T1, class T2, class T = std::common_type_t<T1, T2>, std::enable_if_t<!is_complex_v<T1> && !is_complex_v<T2>, int> = 0>
auto EllipticSNCNDN(std::complex<T1> x, T2 k) -> std::tuple<T, T, T> {
	return EllipticSNCNDN(std::complex<T>(x), T(k));
}


template <class T1, class T2, class = decltype(EllipticSNCNDN(std::declval<const T1&>(), std::declval<const T2&>()))>
auto EllipticSN(T1 x, T2 k) {
	const auto [sn, _1, _2] = EllipticSNCNDN(x, k);
	return sn;
}

template <class T1, class T2, class = decltype(EllipticSNCNDN(std::declval<const T1&>(), std::declval<const T2&>()))>
auto EllipticCN(T1 x, T2 k) {
	const auto [_1, cn, _2] = EllipticSNCNDN(x, k);
	return cn;
}

template <class T1, class T2, class = decltype(EllipticSNCNDN(std::declval<const T1&>(), std::declval<const T2&>()))>
auto EllipticDN(T1 x, T2 k) {
	const auto [_1, _2, dn] = EllipticSNCNDN(x, k);
	return dn;
}

//------------------------------------------------------------------------------
// Jacobi inverse elliptic functions
//------------------------------------------------------------------------------

template <class T, std::enable_if_t<!is_complex_v<T>, int> = 0>
T EllipticArcSN(T x, T k) {
	const auto scale = x * x;
	const auto a = T(1) - x * x;
	const auto b = T(1) - k * k * x * x;
	const auto c = T(1);
	return std::sqrt(scale) * CarlsonRF(a, b, c);
}

template <class T, std::enable_if_t<!is_complex_v<T>, int> = 0>
T EllipticArcCN(T x, T k) {
	const auto scale = T(1) - x * x;
	const auto a = x * x;
	const auto b = T(1);
	const auto c = T(1) - k * k + k * k * x * x;
	return std::sqrt(scale) * CarlsonRF(a, b, c);
}

template <class T, std::enable_if_t<!is_complex_v<T>, int> = 0>
T EllipticArcDN(T x, T k) {
	const auto scale = T(1) - x * x;
	const auto a = k * k * x * x;
	const auto b = k * k;
	const auto c = k * k + x * x - T(1);
	return std::sqrt(scale) * CarlsonRF(a, b, c);
}

template <class T>
std::complex<T> EllipticArcSN(std::complex<T> x, T k) {
	const auto scale = x * x;
	const auto a = T(1) - x * x;
	const auto b = T(1) - k * k * x * x;
	const auto c = T(1);
	return std::sqrt(scale) * CarlsonRF(a, b, c);
}

template <class T>
std::complex<T> EllipticArcCN(std::complex<T> x, T k) {
	const auto scale = T(1) - x * x;
	const auto a = x * x;
	const auto b = T(1);
	const auto c = T(1) - k * k + k * k * x * x;
	return std::sqrt(scale) * CarlsonRF(a, b, c);
}

template <class T>
std::complex<T> EllipticArcDN(std::complex<T> x, T k) {
	const auto scale = T(1) - x * x;
	const auto a = k * k * x * x;
	const auto b = k * k;
	const auto c = k * k + x * x - T(1);
	return std::sqrt(scale) * CarlsonRF(a, b, c);
}

template <class T1, class T2, std::enable_if_t<!is_complex_v<T1> && !is_complex_v<T2>, int> = 0>
auto EllipticArcSN(T1 x, T2 k) {
	using T = std::common_type_t<T1, T2>;
	return EllipticArcSN(T(x), T(k));
}

template <class T1, class T2, std::enable_if_t<!is_complex_v<T1> && !is_complex_v<T2>, int> = 0>
auto EllipticArcCN(T1 x, T2 k) {
	using T = std::common_type_t<T1, T2>;
	return EllipticArcCN(T(x), T(k));
}

template <class T1, class T2, std::enable_if_t<!is_complex_v<T1> && !is_complex_v<T2>, int> = 0>
auto EllipticArcDN(T1 x, T2 k) {
	using T = std::common_type_t<T1, T2>;
	return EllipticArcDN(T(x), T(k));
}

template <class T1, class T2, std::enable_if_t<!is_complex_v<T2>, int> = 0>
auto EllipticArcSN(std::complex<T1> x, T2 k) {
	using T = std::common_type_t<T1, T2>;
	return EllipticArcSN(std::complex<T>(x), T(k));
}

template <class T1, class T2, std::enable_if_t<!is_complex_v<T2>, int> = 0>
auto EllipticArcCN(std::complex<T1> x, T2 k) {
	using T = std::common_type_t<T1, T2>;
	return EllipticArcCN(std::complex<T>(x), T(k));
}

template <class T1, class T2, std::enable_if_t<!is_complex_v<T2>, int> = 0>
auto EllipticArcDN(std::complex<T1> x, T2 k) {
	using T = std::common_type_t<T1, T2>;
	return EllipticArcDN<T>(std::complex<T>(x), T(k));
}


} // namespace dspbb