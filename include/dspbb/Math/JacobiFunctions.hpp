#pragma once

#include <dspbb/Utility/Numbers.hpp>

#include <complex>

namespace dspbb {

// Not sure is these are needed. Probably not.
static_assert(std::numeric_limits<float>::is_iec559);
static_assert(std::numeric_limits<float>::radix == 2);
static_assert(std::numeric_limits<double>::is_iec559);
static_assert(std::numeric_limits<double>::radix == 2);
static_assert(std::numeric_limits<long double>::is_iec559);
static_assert(std::numeric_limits<long double>::radix == 2);


//------------------------------------------------------------------------------
// Lattice transforms
//------------------------------------------------------------------------------

template <class T>
struct LatticeTransform {
	int variant;
	std::complex<T> z;
	std::complex<T> tau;
	std::complex<T> multiplier = T(1);
	std::complex<T> exponent = T(0);
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
		return std::exp(-T(count % 8) * i_v<T> * pi_v<T> / T(4));
	}
	return std::complex<T>(T(1));
}

template <class T>
LatticeTransform<T> ShiftTau(int variant, const std::complex<T>& z, const std::complex<T>& tau) {
	const auto [remainder, count] = ShiftScalar(real(tau));
	const int newVariant = ShiftVariant(variant, count);
	const std::complex<T> multiplier = ShiftMultiplier<T>(variant, count);
	const std::complex<T> newTau = { remainder, imag(tau) };
	return { newVariant, z, newTau, multiplier, T(0) };
}


int InvertVariant(int variant) {
	constexpr std::array<int, 4> permutations = { 1, 4, 3, 2 };
	return permutations[variant - 1];
}

template <class T>
std::pair<std::complex<T>, std::complex<T>> InvertMultiplier(int variant, const std::complex<T>& z, const std::complex<T>& tau) {
	const std::complex<T> factor = std::sqrt(i_v<T> / tau);
	const std::complex<T> exponent = -i_v<T> / tau * z * z / pi_v<T>;

	if (variant == 1) {
		return { -i_v<T> * factor, exponent };
	}
	return { factor, exponent };
}

template <class T>
LatticeTransform<T> InvertTau(int variant, const std::complex<T>& z, const std::complex<T>& tau) {
	const std::complex<T> newTau = -T(1) / tau;
	const std::complex<T> newZ = z * newTau;
	const auto [multiplier, exponent] = InvertMultiplier(variant, z, tau);
	const int newVariant = InvertVariant(variant);
	return { newVariant, newZ, newTau, multiplier, exponent };
}

template <class T>
LatticeTransform<T> RotateTau(int variant, const std::complex<T>& z, const std::complex<T>& tau) {
	const auto [shiftedVariant, shiftedZ, shiftedTau, shiftMultiplier, shiftExponent] = ShiftTau(variant, z, tau);
	const auto [invertedVariant, invertedZ, invertedTau, inverseMultiplier, inverseExponent] = InvertTau(shiftedVariant, shiftedZ, shiftedTau);
	return { invertedVariant, invertedZ, invertedTau, shiftMultiplier * inverseMultiplier, inverseExponent + shiftExponent };
}

//------------------------------------------------------------------------------
// Fourier series approximations
//------------------------------------------------------------------------------

template <class T>
std::complex<T> ThetaSeriesElement(int variant, int n, const std::complex<T>& z, const std::complex<T>& tau, const std::complex<T>& bias) {
	const T sign = T(1 - 2 * (n % 2) * int(variant == 1 || variant == 4));
	const T lambda = T(n) + T(0.5) * T(variant == 1 || variant == 2);
	const auto phase = variant == 1 ? i_v<T> : T(1);
	const auto polarity = variant == 1 ? -T(1) : T(1);

	const auto qExp = lambda * lambda * pi_v<T> * i_v<T> * tau;
	const auto trigExp = i_v<T> * T(2) * lambda * z;

	const auto v1 = std::exp(qExp - trigExp + bias);
	const auto v2 = std::exp(qExp + trigExp + bias);

	return phase * sign * (v1 + polarity * v2);
}

template <class T>
std::complex<T> ThetaSeries(int variant, const std::complex<T>& z, const std::complex<T>& tau, const std::complex<T>& exponent, int iterations) {
	const bool symmetric = variant == 3 || variant == 4;

	const int nFirst = int(symmetric);
	const auto nLast = iterations;

	std::complex<T> accumulator = std::exp(exponent) * T(symmetric);
	for (int n = nLast - 1; n >= nFirst; --n) {
		accumulator += ThetaSeriesElement(variant, n, z, tau, exponent);
	}

	return accumulator;
}

template <class T>
std::complex<T> ContractZ(const std::complex<T>& z) {
	const T re = std::remainder(real(z), T(2) * pi_v<T>);
	return { re, imag(z) };
}

template <class T>
LatticeTransform<T> ReformulateSeries(int variant, const std::complex<T>& z, const std::complex<T>& tau, T threshold, int iterations) {
	LatticeTransform<T> result{ variant, z, tau };
	while (imag(result.tau) < threshold && iterations > 0) {
		const auto transform = RotateTau(result.variant, result.z, result.tau);
		result = { transform.variant,
				   ContractZ(transform.z),
				   transform.tau,
				   transform.multiplier * result.multiplier,
				   transform.exponent + result.exponent };
		--iterations;
	}
	return result;
}

template <class T>
int EstimateSeriesIterations(const std::complex<T>& tau) {
	const T q = std::exp(-pi_v<T> * imag(tau));
	constexpr T epsilon = std::numeric_limits<T>::epsilon();
	return int(T(1) + std::log(epsilon) / std::log(std::max(q, epsilon)));
}

template <class T>
std::complex<T> Theta(int variant, const std::complex<T>& z, const std::complex<T>& tau, int transformIterations = 6, int seriesIterations = 0) {
	assert(1 <= variant && variant <= 4);
	constexpr T threshold = T(0.65); // |exp(i*pi*tau)| = |q| < 0.208
	const LatticeTransform<T> transform = ReformulateSeries(variant, ContractZ(z), tau, threshold, transformIterations);	
	seriesIterations = seriesIterations == 0 ? EstimateSeriesIterations(transform.tau) : seriesIterations;
	return transform.multiplier * ThetaSeries(transform.variant, transform.z, transform.tau, transform.exponent, seriesIterations);
}


//------------------------------------------------------------------------------
// Interface functions
//------------------------------------------------------------------------------

template <class T>
std::complex<T> theta(int variant, std::complex<T> z, std::complex<T> tau) {
	return Theta(variant, z, tau);
}

template <class T>
std::complex<T> theta_1(const std::complex<T>& z, const std::complex<T>& tau) {
	return theta(1, z, tau);
}

template <class T>
std::complex<T> theta_2(const std::complex<T>& z, const std::complex<T>& tau) {
	return theta(2, z, tau);
}

template <class T>
std::complex<T> theta_3(const std::complex<T>& z, const std::complex<T>& tau) {
	return theta(3, z, tau);
}

template <class T>
std::complex<T> theta_4(const std::complex<T>& z, const std::complex<T>& tau) {
	return theta(4, z, tau);
}

} // namespace dspbb