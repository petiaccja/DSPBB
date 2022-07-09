#pragma once

#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "../Utility/TypeTraits.hpp"

#include <complex>
#include <initializer_list>
#include <numeric>
#include <stdexcept>

namespace dspbb {

template <class T>
class Polynomial {
	static_assert(!is_complex_v<T>, "Polynomials can have real coefficients only.");

public:
	Polynomial() noexcept = default;
	Polynomial(std::initializer_list<T> coefficients);
	template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter>()), T>, int> = 0>
	Polynomial(Iter firstCoefficient, Iter lastCoefficient);

	void resize(size_t numCoefficients, T value = T{});
	size_t size() const;
	size_t order() const;

	BasicSignalView<const T, DOMAINLESS> coefficients() const;
	BasicSignalView<T, DOMAINLESS> coefficients();

	T operator()(const T& x) const;
	std::complex<T> operator()(const std::complex<T>& x) const;

private:
	template <class U>
	U Eval(const U& x) const;

private:
	BasicSignal<T, DOMAINLESS> m_coefficients;
};


template <class T>
class FactoredPolynomial {
	static_assert(!is_complex_v<T>, "Factored polynomials are templated with real numbers only.");

public:
	FactoredPolynomial() noexcept = default;
	FactoredPolynomial(FactoredPolynomial&& rhs) noexcept;
	FactoredPolynomial(const FactoredPolynomial& rhs);
	FactoredPolynomial& operator=(FactoredPolynomial&& rhs) noexcept;
	FactoredPolynomial& operator=(const FactoredPolynomial& rhs);
	~FactoredPolynomial() = default;

	FactoredPolynomial(std::initializer_list<std::complex<T>> roots);
	template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter>()), std::complex<T>>, int> = 0>
	FactoredPolynomial(Iter firstRoot, Iter lastRoot);

	void resize(size_t numRealRoots, size_t numComplexPairs, T realValue = T{}, std::complex<T> complexValue = std::complex<T>{});
	void regroup(size_t numRealRoots, T realValue = T{}, std::complex<T> complexValue = std::complex<T>{});
	size_t num_roots() const { return num_real_roots() + num_complex_roots(); }
	size_t num_real_roots() const { return m_real.size(); }
	size_t num_complex_roots() const { return 2 * num_complex_pairs(); }
	size_t num_complex_pairs() const { return m_complex.size(); }
	size_t order() const { return num_roots(); }

	BasicSignalView<const T, DOMAINLESS> real_roots() const;
	BasicSignalView<const std::complex<T>, DOMAINLESS> complex_pairs() const;
	BasicSignalView<T, DOMAINLESS> real_roots();
	BasicSignalView<std::complex<T>, DOMAINLESS> complex_pairs();

	T operator()(const T& x) const;
	std::complex<T> operator()(const std::complex<T>& x) const;

private:
	template <class U>
	U eval(const U& x) const;

private:
	BasicSignal<T, DOMAINLESS> m_mem;
	BasicSignalView<T, DOMAINLESS> m_real;
	BasicSignalView<std::complex<T>, DOMAINLESS> m_complex;
};


template <class T>
Polynomial<T>::Polynomial(std::initializer_list<T> coefficients) : m_coefficients(coefficients) {}

template <class T>
template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter>()), T>, int>>
Polynomial<T>::Polynomial(Iter firstCoefficient, Iter lastCoefficient) : m_coefficients(firstCoefficient, lastCoefficient) {}

template <class T>
void Polynomial<T>::resize(size_t numCoefficients, T value) {
	m_coefficients.resize(numCoefficients, value);
}

template <class T>
size_t Polynomial<T>::size() const {
	return m_coefficients.size();
}

template <class T>
size_t Polynomial<T>::order() const {
	return size() != 0 ? size() - 1 : 0;
}

template <class T>
BasicSignalView<const T, DOMAINLESS> Polynomial<T>::coefficients() const {
	return AsView(m_coefficients);
}

template <class T>
BasicSignalView<T, DOMAINLESS> Polynomial<T>::coefficients() {
	return AsView(m_coefficients);
}

template <class T>
T Polynomial<T>::operator()(const T& x) const {
	return Eval(x);
}

template <class T>
std::complex<T> Polynomial<T>::operator()(const std::complex<T>& x) const {
	return Eval(x);
}

template <class T>
template <class U>
U Polynomial<T>::Eval(const U& x) const {
	U y = U(0);
	for (auto it = m_coefficients.rbegin(); it != m_coefficients.rend(); ++it) {
		const auto& c = *it;
		y = x * y + c;
	}
	return y;
}


namespace impl {

	template <class Iter>
	bool AreRootsConjugatePairs(Iter first, Iter last) {
		return std::all_of(first, last, [&first, &last](const auto& root) {
			size_t numSame = std::count_if(first, last, [&root](const auto& other) { return other == root; });
			size_t numConj = std::count_if(first, last, [&root](const auto& other) { return other == std::conj(root); });
			return numSame == numConj;
		});
	}

} // namespace impl

template <class T>
FactoredPolynomial<T>::FactoredPolynomial(FactoredPolynomial&& rhs) noexcept
	: m_mem{ std::move(rhs.m_mem) },
	  m_real{ m_mem.data(), rhs.m_real.size() },
	  m_complex{ reinterpret_cast<std::complex<T>*>(m_mem.data() + rhs.m_real.size()), rhs.m_complex.size() } {
	rhs.m_real = {};
	rhs.m_complex = {};
}


template <class T>
FactoredPolynomial<T>::FactoredPolynomial(const FactoredPolynomial& rhs)
	: m_mem{ rhs.m_mem },
	  m_real{ m_mem.data(), rhs.m_real.size() },
	  m_complex{ reinterpret_cast<std::complex<T>*>(m_mem.data() + rhs.m_real.size()), rhs.m_complex.size() } {}

template <class T>
FactoredPolynomial<T>& FactoredPolynomial<T>::operator=(FactoredPolynomial&& rhs) noexcept {
	m_mem = std::move(rhs.m_mem),
	m_real = { m_mem.data(), rhs.m_real.size() },
	m_complex = { reinterpret_cast<std::complex<T>*>(m_mem.data() + rhs.m_real.size()), rhs.m_complex.size() };
	rhs.m_real = {};
	rhs.m_complex = {};
}

template <class T>
FactoredPolynomial<T>& FactoredPolynomial<T>::operator=(const FactoredPolynomial& rhs) {
	m_mem = rhs.m_mem,
	m_real = { m_mem.data(), rhs.m_real.size() },
	m_complex = { reinterpret_cast<std::complex<T>*>(m_mem.data() + rhs.m_real.size()), rhs.m_complex.size() };
}


template <class T>
FactoredPolynomial<T>::FactoredPolynomial(std::initializer_list<std::complex<T>> roots) : FactoredPolynomial(roots.begin(), roots.end()) {}


template <class T>
template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter>()), std::complex<T>>, int>>
FactoredPolynomial<T>::FactoredPolynomial(Iter firstRoot, Iter lastRoot) {
	const size_t numReal = std::transform_reduce(firstRoot, lastRoot, size_t(0), std::plus<size_t>{}, [](const auto& arg) { return size_t(imag(arg) == T(0)); });
	const size_t numComplex = std::distance(firstRoot, lastRoot) - numReal;
	if (numComplex % 2 != 0) {
		throw std::invalid_argument{ "Both of the complex conjugate pair roots have to be provided." };
	}
	if (!impl::AreRootsConjugatePairs(firstRoot, lastRoot)) {
		throw std::invalid_argument{ "All complex roots must form conjugate pairs." };
	}

	resize(numReal, numComplex / 2);
	auto realIt = m_real.begin();
	for (auto rootIt = firstRoot; rootIt != lastRoot; ++rootIt) {
		if (imag(*rootIt) == T(0)) {
			*realIt = real(*rootIt);
			++realIt;
		}
	}
	std::copy_if(firstRoot, lastRoot, m_complex.begin(), [](const auto& arg) { return imag(arg) > T(0); });
}

template <class T>
void FactoredPolynomial<T>::resize(size_t numRealRoots, size_t numComplexPairs, T realValue, std::complex<T> complexValue) {
	const size_t numCoeffs = numRealRoots + 2 * numComplexPairs;
	const size_t currentRealRoots = num_real_roots();
	const size_t currentComplexPairs = num_complex_pairs();
	if (numRealRoots < currentRealRoots) {
		std::move(m_mem.begin() + currentRealRoots, m_mem.end(), m_mem.begin() + numRealRoots);
	}
	m_mem.resize(numCoeffs);
	if (currentRealRoots < numRealRoots) {
		std::move_backward(m_mem.begin() + currentRealRoots, m_mem.begin() + (numCoeffs - numRealRoots), m_mem.end());
	}
	m_real = { m_mem.data(), numRealRoots };
	m_complex = { reinterpret_cast<std::complex<T>*>(m_mem.data() + numRealRoots), numComplexPairs };
	if (currentRealRoots < numRealRoots) {
		std::fill(m_real.begin() + currentRealRoots, m_real.end(), realValue);
	}
	if (currentComplexPairs < numComplexPairs) {
		std::fill(m_complex.begin() + currentComplexPairs, m_complex.end(), complexValue);
	}
}

template <class T>
void FactoredPolynomial<T>::regroup(size_t numRealRoots, T realValue, std::complex<T> complexValue) {
	if (numRealRoots % 2 != num_roots() % 2) {
		throw std::invalid_argument("You can't have complex roots that are not conjugate pairs.");
	}
	if (numRealRoots > num_roots()) {
		throw std::invalid_argument("Regrouping must preserve the total number of roots.");
	}

	const size_t numComplexPairs = (num_roots() - numRealRoots) / 2;
	const size_t currentRealRoots = num_real_roots();
	const size_t currentComplexPairs = num_complex_pairs();
	if (numRealRoots < currentRealRoots) {
		std::move(m_mem.begin() + currentRealRoots, m_mem.end(), m_mem.begin() + numRealRoots);
	}
	if (currentRealRoots < numRealRoots) {
		std::move_backward(m_mem.begin() + currentRealRoots, m_mem.end() - (numRealRoots - currentRealRoots), m_mem.end());
	}
	m_real = { m_mem.data(), numRealRoots };
	m_complex = { reinterpret_cast<std::complex<T>*>(m_mem.data() + numRealRoots), numComplexPairs };
	if (currentRealRoots < numRealRoots) {
		std::fill(m_real.begin() + currentRealRoots, m_real.end(), realValue);
	}
	if (currentComplexPairs < numComplexPairs) {
		std::fill(m_complex.begin() + currentComplexPairs, m_complex.end(), complexValue);
	}
}

template <class T>
BasicSignalView<const T, DOMAINLESS> FactoredPolynomial<T>::real_roots() const {
	return AsView(m_real);
}

template <class T>
BasicSignalView<const std::complex<T>, DOMAINLESS> FactoredPolynomial<T>::complex_pairs() const {
	return AsView(m_complex);
}

template <class T>
BasicSignalView<T, DOMAINLESS> FactoredPolynomial<T>::real_roots() {
	return AsView(m_real);
}

template <class T>
BasicSignalView<std::complex<T>, DOMAINLESS> FactoredPolynomial<T>::complex_pairs() {
	return AsView(m_complex);
}

template <class T>
T FactoredPolynomial<T>::operator()(const T& x) const {
	return eval(x);
}

template <class T>
std::complex<T> FactoredPolynomial<T>::operator()(const std::complex<T>& x) const {
	return eval(x);
}

template <class T>
template <class U>
U FactoredPolynomial<T>::eval(const U& x) const {
	const U rp = std::transform_reduce(m_real.begin(), m_real.end(), U(T(1)), std::multiplies{}, [&x](const auto& root) {
		return x - root;
	});
	const U cp = std::transform_reduce(m_complex.begin(), m_complex.end(), U(T(1)), std::multiplies{}, [&x](const auto& root) {
		const auto a = real(root);
		const auto b = imag(root);
		return x * x - T(2) * a * x + a * a + b * b;
	});
	return rp * cp;
}


namespace impl {
	template <class T>
	void MultiplyPolynomialBy1stOrder(BasicSignalView<std::complex<T>, DOMAINLESS> coefficients, std::complex<T> c0) {
		for (auto it1 = coefficients.rbegin(); it1 != coefficients.rend(); ++it1) {
			auto it0 = it1;
			++it0;
			const std::complex<T> p0 = it0 != coefficients.rend() ? *it0 : 0.0f;
			const std::complex<T> p1 = *it1;
			*it1 = p0 + p1 * c0;
		}
	}

	template <class T>
	void MultiplyPolynomialBy1stOrder(BasicSignalView<T, DOMAINLESS> coefficients, T c0) {
		for (auto it1 = coefficients.rbegin(); it1 != coefficients.rend(); ++it1) {
			auto it0 = it1;
			++it0;
			const T p0 = it0 != coefficients.rend() ? *it0 : 0.0f;
			const T p1 = *it1;
			*it1 = p0 + p1 * c0;
		}
	}

	template <class T>
	void MultiplyPolynomialBy2ndOrder(BasicSignalView<T, DOMAINLESS> coefficients, T c0, T c1) {
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
Polynomial<T> ExpandPolynomial(const FactoredPolynomial<T>& factored) {
	Polynomial<T> poly;
	poly.resize(factored.num_roots() + 1, T(0));
	poly.coefficients()[0] = T(1);
	for (const auto& root : factored.real_roots()) {
		impl::MultiplyPolynomialBy1stOrder(poly.coefficients(), -root);
	}
	for (const auto& root : factored.complex_pairs()) {
		const T real = std::real(root);
		const T imag = std::imag(root);
		impl::MultiplyPolynomialBy2ndOrder(poly.coefficients(), real * real + imag * imag, -2.0f * real);
	}
	return poly;
}


} // namespace dspbb