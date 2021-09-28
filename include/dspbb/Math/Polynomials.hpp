#pragma once

#include <complex>
#include <vector>

namespace dspbb {

namespace impl {
	template <class T>
	void MultiplyPolynomialBy1stOrder(std::vector<std::complex<T>>& coefficients, std::complex<T> c0) {
		for (auto it1 = coefficients.rbegin(); it1 != coefficients.rend(); ++it1) {
			auto it0 = it1;
			++it0;
			const std::complex<T> p0 = it0 != coefficients.rend() ? *it0 : 0.0f;
			const std::complex<T> p1 = *it1;
			*it1 = p0 + p1 * c0;
		}
	}

	template <class T>
	void MultiplyPolynomialBy1stOrder(std::vector<T>& coefficients, T c0) {
		for (auto it1 = coefficients.rbegin(); it1 != coefficients.rend(); ++it1) {
			auto it0 = it1;
			++it0;
			const T p0 = it0 != coefficients.rend() ? *it0 : 0.0f;
			const T p1 = *it1;
			*it1 = p0 + p1 * c0;
		}
	}

	template <class T>
	void MultiplyPolynomialBy2ndOrder(std::vector<T>& coefficients, T c0, T c1) {
		for (auto it2 = coefficients.rbegin(); it2 != coefficients.rend(); ++it2) {
			auto it1 = it2;
			++it1;
			auto it0 = it1;
			if (it0 != coefficients.rend()) {
				++it0;
			}
			const T p0 = it0 != coefficients.rend() ? *it0 : 0.0f;
			const T p1 = it1 != coefficients.rend() ? *it1 : 0.0f;
			const T p2 = *it2;
			*it2 = p0 + p1 * c1 + p2 * c0;
		}
	}
} // namespace impl

template <class T>
std::vector<std::complex<T>> ExpandPolynomialComplex(const std::vector<std::complex<T>>& roots) {
	std::vector<std::complex<T>> coefficients;
	coefficients.resize(roots.size() + 1, 0.0f);
	coefficients[0] = 1.0f;
	for (const auto& root : roots) {
		impl::MultiplyPolynomialBy1stOrder(coefficients, -root);
	}
	return coefficients;
}

template <class T>
std::vector<T> ExpandPolynomialReal(const std::vector<std::complex<T>>& roots) {
	std::vector<T> coefficients;
	coefficients.resize(roots.size() + 1, 0.0f);
	coefficients[0] = 1.0f;
	for (const auto& root : roots) {
		const T real = std::real(root);
		const T imag = std::imag(root);
		if (imag > 0.0f) {
			impl::MultiplyPolynomialBy2ndOrder(coefficients, real * real + imag * imag, -2.0f * real);
		}
		else if (imag == 0.0f) {
			impl::MultiplyPolynomialBy1stOrder(coefficients, -real);
		}
	}
	return coefficients;
}

template <class T, class U>
auto EvaluatePolynomial(const std::vector<T>& coefficients, U x) {
	U xpow = static_cast<U>(1);
	auto acc = T(0) * U(0);
	for (const auto& coeff : coefficients) {
		acc += xpow * coeff;
		xpow *= x;
	}
	return acc;
}


} // namespace dspbb