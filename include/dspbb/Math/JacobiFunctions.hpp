#pragma once

#include <dspbb/Utility/Numbers.hpp>

#include <complex>

namespace dspbb {

static_assert(std::numeric_limits<float>::is_iec559);
static_assert(std::numeric_limits<float>::radix == 2);
static_assert(std::numeric_limits<double>::is_iec559);
static_assert(std::numeric_limits<double>::radix == 2);
static_assert(std::numeric_limits<long double>::is_iec559);
static_assert(std::numeric_limits<long double>::radix == 2);

template <class T>
constexpr std::complex<T> i_v = { T(0), T(1) };


template <class T>
struct TauTransform {
	std::complex<T> z;
	std::complex<T> tau;
	std::complex<T> multiplier;
	std::complex<T> exponent;
	int variant;
};


template <class T>
std::pair<T, int> ShiftScalar(T value) {
	int integer;
	T fractional = std::remquo(value, T(1), &integer);
	if (fractional < -T(0.5)) {
		fractional += T(1);
		integer += 1;
	}
	else if (fractional > T(0.5)) {
		fractional -= T(1);
		integer -= 1;
	}
	return { fractional, -integer };
}

int ShiftVariant(int variant, int count) {
	constexpr std::array<int, 4> permutations = { 1, 2, 4, 3 };
	return count % 2 == 0 ? variant : permutations[variant - 1];
}

template <class T>
std::complex<T> ShiftMultiplier(int variant, int count) {
	if (variant == 1 || variant == 2) {
		return std::exp(std::complex<T>(T(0), T(1)) * pi_v<T> * T(count % 8) / T(4));
	}
	return std::complex<T>(T(1));
}

template <class T>
TauTransform<T> ShiftTau(const std::complex<T>& tau, int variant) {
	const auto [remainder, count] = ShiftScalar(real(tau));
	const int newVariant = ShiftVariant(variant, count);
	const std::complex<T> multiplier = ShiftMultiplier<T>(variant, count);
	const std::complex<T> newTau = { remainder, imag(tau) };
	return { {}, newTau, multiplier, { T(0) }, newVariant };
}


int InvertVariant(int variant) {
	constexpr std::array<int, 4> permutations = { 1, 4, 3, 2 };
	return permutations[variant - 1];
}

template <class T>
std::pair<std::complex<T>, std::complex<T>> InvertMultiplier(const std::complex<T> z, const std::complex<T> tau, int variant) {
	const std::complex<T> den = std::sqrt(-i_v<T> * tau);
	const std::complex<T> exponent = -i_v<T> * T(1) / tau * z * z / pi_v<T>;

	if (variant == 1) {
		return { -i_v<T> / den, exponent };
	}
	else {
		return { T(1) / den, exponent };
	}
}

template <class T>
TauTransform<T> InvertTau(const std::complex<T>& z, const std::complex<T>& tau, int variant) {
	const std::complex<T> newTau = -T(1) / tau;
	const std::complex<T> newZ = z * newTau;
	const auto [multiplier, exponent] = InvertMultiplier(z, tau, variant);
	const int newVariant = InvertVariant(variant);
	return { newZ, newTau, multiplier, exponent, newVariant };
}


template <class T>
TauTransform<T> RotateTau(const std::complex<T>& z, const std::complex<T>& tau, int variant) {
	const auto [_, shiftedTau, shiftMultiplier, __, shiftedVariant] = ShiftTau(tau, variant);
	const auto [invertedZ, invertedTau, inverseMultiplier, inverseExponent, invertedVariant] = InvertTau(z, shiftedTau, shiftedVariant);
	return { invertedZ, invertedTau, shiftMultiplier * inverseMultiplier, inverseExponent, invertedVariant };
}


template <class T>
std::complex<T> ThetaSeriesElement1(int n, const std::complex<T>& z, const std::complex<T>& tau, const std::complex<T>& bias) {
	const T sign = static_cast<T>(1 - 2 * (n % 2));
	const T exponent = static_cast<T>(n) + T(0.5);

	const auto u = exponent * exponent * pi_v<T> * i_v<T> * tau;
	const auto v = i_v<T> * T(2 * n + 1) * z;
	return i_v<T> * sign * (std::exp(u - v + bias) - std::exp(u + v + bias));
}

template <class T>
std::complex<T> ThetaSeriesElement2(int n, const std::complex<T>& z, const std::complex<T>& tau, const std::complex<T>& bias) {
	const T exponent = static_cast<T>(n) + T(0.5);

	const auto u = exponent * exponent * pi_v<T> * i_v<T> * tau;
	const auto v = i_v<T> * T(2 * n + 1) * z;
	return std::exp(u - v + bias) + std::exp(u + v + bias);
}

template <class T>
std::complex<T> ThetaSeriesElement3(int n, const std::complex<T>& z, const std::complex<T>& tau, const std::complex<T>& bias) {
	const T exponent = static_cast<T>(n);

	const auto u = exponent * exponent * pi_v<T> * i_v<T> * tau;
	const auto v = i_v<T> * T(2 * n) * z;
	return std::exp(u - v + bias) + std::exp(u + v + bias);
}

template <class T>
std::complex<T> ThetaSeriesElement4(int n, const std::complex<T>& z, const std::complex<T>& tau, const std::complex<T>& bias) {
	const T sign = static_cast<T>(1 - 2 * (n % 2));
	const T exponent = static_cast<T>(n);

	const auto u = exponent * exponent * pi_v<T> * i_v<T> * tau;
	const auto v = i_v<T> * T(2 * n) * z;
	return sign * (std::exp(u - v + bias) + std::exp(u + v + bias));
}

template <class T>
std::complex<T> ThetaSeries(std::complex<T> z, std::complex<T> tau, const std::complex<T>& exponent, int variant) {
	constexpr auto numIters = 5;

	bool symmetricVariant = variant == 3 || variant == 4;
	std::complex<T> acc{ T(symmetricVariant), T(0) };

	for (int n = int(symmetricVariant); n < numIters; ++n) {
		std::complex<T> element;
		switch (variant) {
			case 1: element = ThetaSeriesElement1(n, z, tau, exponent); break;
			case 2: element = ThetaSeriesElement2(n, z, tau, exponent); break;
			case 3: element = ThetaSeriesElement3(n, z, tau, exponent); break;
			case 4: element = ThetaSeriesElement4(n, z, tau, exponent); break;
		}
		acc += element;
	}

	return acc;
}


template <class T>
std::complex<T> theta(std::complex<T> z, std::complex<T> tau, int variant) {
	constexpr int maxIter = 0;
	int iter = 0;
	std::complex<T> multiplier = { T(1), T(0) };
	std::complex<T> exponent = { T(0), T(0) };
	while (imag(tau) < T(0.5) && iter++ < maxIter) {
		const auto [newZ, newTau, newMultiplier, newExponent, newVariant] = RotateTau(z, tau, variant);
		z = newZ;
		tau = newTau;
		variant = newVariant;
		multiplier *= newMultiplier;
		exponent += newExponent;
	}
	if (iter == maxIter && iter > 0) {
		throw std::domain_error("tau should be on the upper half-plane. Algorithm may fail if Im[tau] is very small.");
	}
	return multiplier * ThetaSeries(z, tau, exponent, variant);
}


template <class T>
std::complex<T> theta_1(const std::complex<T>& z, const std::complex<T>& tau) {
	return theta(z, tau, 1);
}

template <class T>
std::complex<T> theta_2(const std::complex<T>& z, const std::complex<T>& tau) {
	return theta(z, tau, 2);
}

template <class T>
std::complex<T> theta_3(const std::complex<T>& z, const std::complex<T>& tau) {
	return theta(z, tau, 3);
}

template <class T>
std::complex<T> theta_4(const std::complex<T>& z, const std::complex<T>& tau) {
	return theta(z, tau, 4);
}

} // namespace dspbb