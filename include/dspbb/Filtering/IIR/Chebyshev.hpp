#pragma once

#include "../../Generators/Spaces.hpp"
#include "../../LTISystems/Systems.hpp"
#include "../../Utility/Numbers.hpp"


namespace dspbb {

namespace impl {

	template <class T>
	std::complex<T> Chebyshev1Pole(size_t index, size_t order, T epsilon) {
		constexpr std::complex<T> i{ T(0), T(1) };
		const std::complex<T> theta = T(1) / T(order) * std::acos(+i / epsilon) + T(index) / T(order) * pi_v<T>;
		return i * std::cos(theta);
	}

	template <class T>
	std::complex<T> Chebyshev2Pole(size_t index, size_t order, T epsilon) {
		return T(1) / Chebyshev1Pole(index, order, epsilon);
	}

	template <class T>
	std::complex<T> Chebyshev2Zero(size_t index, size_t order, T epsilon) {
		constexpr std::complex<T> i{ T(0), T(1) };
		const std::complex<T> zero = -i * std::cos(pi_v<T> / T(2) * T(2 * index + 1) / T(order));
		return T(1) / zero;
	}


} // namespace impl

template <class T>
ZeroPoleGain<T, eDiscretization::CONTINUOUS> Chebyshev1(size_t order, T ripple) {
	const T epsilon = std::sqrt(T(1) / std::pow(T(1) - ripple, T(2)) - 1);
	const T gain = T(1) / (T(1u << (order - 1)) * epsilon);

	FactoredPolynomial<T> poles;
	poles.Resize(order % 2, order / 2);

	size_t index = 0;
	for (auto& root : poles.ComplexPairs()) {
		root = impl::Chebyshev1Pole(index, order, epsilon);
		++index;
	}
	for (auto& root : poles.RealRoots()) {
		root = std::real(impl::Chebyshev1Pole(index, order, epsilon));
	}

	return { gain, {}, std::move(poles) };
}


template <class T>
ZeroPoleGain<T, eDiscretization::CONTINUOUS> Chebyshev2(size_t order, T ripple) {
	const T epsilon = ripple / std::sqrt(T(1) - ripple * ripple);
	const T gain = order % 2 == 0 ? ripple : epsilon * T(order);

	FactoredPolynomial<T> poles;
	FactoredPolynomial<T> zeros;
	poles.Resize(order % 2, order / 2);
	zeros.Resize(0, order / 2);

	std::vector<std::complex<T>> p;
	std::vector<std::complex<T>> z;

	size_t index = 0;
	for (auto& root : poles.ComplexPairs()) {
		root = impl::Chebyshev2Pole(index, order, epsilon);
		p.push_back(root);
		++index;
	}
	for (auto& root : poles.RealRoots()) {
		root = std::real(impl::Chebyshev2Pole(index, order, epsilon));
		p.push_back(root);
	}

	index = 0;
	for (auto& root : zeros.ComplexPairs()) {
		root = impl::Chebyshev2Zero(index, order, epsilon);
		z.push_back(root);
		++index;
	}

	return { gain, std::move(zeros), std::move(poles) };
}


} // namespace dspbb