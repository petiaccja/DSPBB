#pragma once

#include <complex>
#include <vector>
#include <initializer_list>

namespace dspbb {

template <class T>
class Polynomial {
	static_assert(!is_complex_v<T>, "Polynomials can have real coefficients only.");

public:
	Polynomial() = default;
	Polynomial(std::initializer_list<T> coefficients);
	template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter>()), T>, int> = 0>
	Polynomial(Iter firstCoefficient, Iter lastCoefficient);

	SignalView<const T, DOMAINLESS> Coefficients() const;
	SignalView<T, DOMAINLESS> Coefficients();

	T operator()(const T& x) const;
	std::complex<T> operator()(const std::complex<T>& x) const;

private:
	Signal<T, DOMAINLESS> m_coefficients;
};


template <class T>
class FactoredPolynomial {
	static_assert(!is_complex_v<T>, "Factored polynomials are templated with real numbers only.");

public:
	FactoredPolynomial() = default;
	FactoredPolynomial(std::initializer_list<std::complex<T>> roots);
	template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter>()), std::complex<T>>, int> = 0>
	FactoredPolynomial(Iter firstRoot, Iter lastRoot);

	SignalView<const T, DOMAINLESS> RealRoots() const;
	SignalView<const std::complex<T>, DOMAINLESS> ComplexRoots() const;
	SignalView<T, DOMAINLESS> RealRoots();
	SignalView<std::complex<T>, DOMAINLESS> ComplexRoots();

	T operator()(const T& x) const;
	std::complex<T> operator()(const std::complex<T>& x) const;

private:
	Signal<T, DOMAINLESS> m_mem;
	SignalView<T, DOMAINLESS> m_real;
	SignalView<std::complex<T>, DOMAINLESS> m_complex;
};



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