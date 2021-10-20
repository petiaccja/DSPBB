#pragma once

#include "../Primitives/Signal.hpp"
#include "../Utility/TypeTraits.hpp"
#include "System.hpp"

#include <numeric>
#include <vector>


namespace dspbb {


template <class T, eSystemDiscretization Discretization>
struct PoleZeroSystem {
	static_assert(!is_complex_v<T>, "Systems are templated with real numbers only.");

public:
	using C = std::complex<T>;
	PoleZeroSystem(T gain = 1, std::vector<C> poles = {}, std::vector<C> zeros = {});
	const std::vector<C>& Poles() const { return m_poles; }
	const std::vector<C>& Zeros() const { return m_zeros; }
	const T& Gain() const { return m_gain; }
	C operator()(const C& x) const;
	T operator()(const T& x) const;

private:
	T m_gain;
	std::vector<C> m_poles;
	std::vector<C> m_zeros;
};

template <class T>
using ContinuousPoleZeroSystem = PoleZeroSystem<T, eSystemDiscretization::CONTINUOUS>;

template <class T>
using DiscretePoleZeroSystem = PoleZeroSystem<T, eSystemDiscretization::DISCRETE>;


namespace impl {

	template <class Iter>
	bool CheckConjugateRoots(Iter first, Iter last) {
		return std::all_of(first, last, [&first, &last](const auto& root) {
			size_t numSame = std::count_if(first, last, [&root](const auto& other) { return other == root; });
			size_t numConj = std::count_if(first, last, [&root](const auto& other) { return other == std::conj(root); });
			return numSame == numConj;
		});
	}

} // namespace impl


template <class T, eSystemDiscretization Discretization>
PoleZeroSystem<T, Discretization>::PoleZeroSystem(T gain, std::vector<C> poles, std::vector<C> zeros)
	: m_gain(gain),
	  m_poles(std::move(poles)),
	  m_zeros(std::move(zeros)) {
	if (!impl::CheckConjugateRoots(m_poles.begin(), m_poles.end()) || !impl::CheckConjugateRoots(m_zeros.begin(), m_zeros.end())) {
		throw std::invalid_argument("Every complex pole must be an exact (floating point) conjugate pair.");
	}
}

template <class T, eSystemDiscretization Discretization>
auto PoleZeroSystem<T, Discretization>::operator()(const C& x) const -> C {
	const auto EvalFactoredPolynomial = [](const auto& x, auto first, auto last) {
		return std::accumulate(first, last, C(1), [&x](const C& acc, const C& item) { return acc * (x - item); });
	};
	const auto num = EvalFactoredPolynomial(x, m_zeros.begin(), m_zeros.end());
	const auto den = EvalFactoredPolynomial(x, m_poles.begin(), m_poles.end());
	return m_gain * num / den;
}

template <class T, eSystemDiscretization Discretization>
T PoleZeroSystem<T, Discretization>::operator()(const T& x) const {
	const auto complexResult = operator()(C{ x }); // Imaginary part should be zero due to conjugate symmetry.
	return std::real(complexResult);
}


} // namespace dspbb
