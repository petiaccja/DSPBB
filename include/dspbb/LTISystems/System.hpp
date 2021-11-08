#pragma once

#include "../Math/Polynomials.hpp"


namespace dspbb {

enum class eDiscretization {
	DISCRETE,
	CONTINUOUS,
};


template <class T, eDiscretization Discretization>
struct ZeroPoleGain {
	T gain;
	FactoredPolynomial<T> zeros;
	FactoredPolynomial<T> poles;

	std::complex<T> operator()(const std::complex<T>& x) const { return gain * zeros(x) / poles(x); }
	T operator()(const T& x) const { return gain * zeros(x) / poles(x); }
};


template <class T, eDiscretization Discretization>
struct TransferFunction {
	TransferFunction() = default;
	TransferFunction(const ZeroPoleGain<T, Discretization>& zpk);

	Polynomial<T> numerator;
	Polynomial<T> denominator;

	std::complex<T> operator()(const std::complex<T>& x) const { return numerator(x) / denominator(x); }
	T operator()(const T& x) const { return numerator(x) / denominator(x); }
};


template <class T>
struct CascadedBiquad {
	CascadedBiquad() = default;
	CascadedBiquad(const ZeroPoleGain<T, eDiscretization::DISCRETE>& zpk);

	struct Biquad {
		std::array<T, 3> numerator = { 1, 0, 0 };
		std::array<T, 2> denominator = { 1, 0 };
	};
	std::vector<Biquad> sections;

private:
	template <class Iter>
	std::vector<std::array<T, 3>> RealRootPolynomials(Iter first, Iter last);
	template <class Iter>
	std::vector<std::array<T, 3>> ComplexPairPolynomials(Iter first, Iter last);
};


template <class T, eDiscretization Discretization>
TransferFunction<T, Discretization>::TransferFunction(const ZeroPoleGain<T, Discretization>& zpk)
	: numerator{ ExpandPolynomial(zpk.zeros) },
	  denominator{ ExpandPolynomial(zpk.poles) } {
	numerator.Coefficients() *= zpk.gain;
}

template <class T>
CascadedBiquad<T>::CascadedBiquad(const ZeroPoleGain<T, eDiscretization::DISCRETE>& zpk) {
	const auto realZeroPolys = RealRootPolynomials(zpk.zeros.RealRoots().begin(), zpk.zeros.RealRoots().end());
	const auto complexZeroPolys = ComplexPairPolynomials(zpk.zeros.ComplexPairs().begin(), zpk.zeros.ComplexPairs().end());
	const auto realPolePolys = RealRootPolynomials(zpk.poles.RealRoots().begin(), zpk.poles.RealRoots().end());
	const auto complexPolePolys = ComplexPairPolynomials(zpk.poles.ComplexPairs().begin(), zpk.poles.ComplexPairs().end());

	auto zeroPolys = realZeroPolys;
	zeroPolys.insert(complexZeroPolys.begin(), complexZeroPolys.end());
	auto polePolys = realPolePolys;
	polePolys.insert(realPolePolys.begin(), realPolePolys.end());

	const auto polyAscending = [](const auto& lhs, const auto& rhs) {
		return std::abs(lhs[0]) < std::abs(rhs[0]);
	};
	std::sort(zeroPolys.begin(), zeroPolys.end(), polyAscending);
	std::sort(polePolys.begin(), polePolys.end(), polyAscending);

	const size_t numSections = std::max(zeroPolys.size(), polePolys.size());
	sections.resize(numSections);
	for (size_t i = 0; i < zeroPolys.size(); ++i) {
		sections[i].numerator = zeroPolys[i];
	}
	for (size_t i = 0; i < polePolys.size(); ++i) {
		sections[i].denominator = { polePolys[i][0], polePolys[i][1] };
	}
}

template <class T>
template <class Iter>
std::vector<std::array<T, 3>> CascadedBiquad<T>::RealRootPolynomials(Iter first, Iter last) {
	std::vector<T> ascending{ first, last };
	std::sort(ascending.begin(), ascending.end());

	std::vector<std::array<T, 3>> polynomials;

	auto firstPositive = std::partition_point(ascending.begin(), ascending.end(), [](const T& r) { r < T(0); });
	auto lastPositive = ascending.end();
	auto firstNegative = ascending.begin();
	auto lastNegative = firstPositive != ascending.begin() ? firstPositive - 1 : firstPositive;
	for (; firstNegative < lastNegative; ++firstNegative, --lastNegative) {
		const auto& r1 = *firstNegative;
		const auto& r2 = *lastNegative;
		polynomials.push_back({ r1 * r2, -r1 - r2, T(1) });
	}
	for (; firstPositive < lastPositive; ++firstPositive, --lastPositive) {
		const auto& r1 = *firstNegative;
		const auto& r2 = *lastNegative;
		polynomials.push_back({ r1 * r2, -r1 - r2, T(1) });
	}
	const bool unpairedPositive = firstPositive == lastPositive;
	const bool unpairedNegative = firstNegative == lastNegative;
	if (unpairedNegative && unpairedPositive) {
		const auto& r1 = *firstPositive;
		const auto& r2 = *firstNegative;
		polynomials.push_back({ r1 * r2, -r1 - r2, T(1) });
	}
	else if (unpairedPositive) {
		polynomials.push_back({ *firstPositive, T(1), T(0) });
	}
	else if (unpairedNegative) {
		polynomials.push_back({ *firstNegative, T(1), T(0) });
	}

	return polynomials;
}

template <class T>
template <class Iter>
std::vector<std::array<T, 3>> CascadedBiquad<T>::ComplexPairPolynomials(Iter first, Iter last) {
	std::vector<std::array<T, 3>> polynomials;
	while (first != last) {
		std::complex<T> r = *first;
		polynomials.push_back({ std::abs(r) * std::abs(r), -T(2) * std::real(r), T(1) });
	}
	return polynomials;
}


template <class T>
using ContinuousTransferFunction = TransferFunction<T, eDiscretization::CONTINUOUS>;
template <class T>
using DiscreteTransferFunction = TransferFunction<T, eDiscretization::DISCRETE>;

template <class T>
using ContinuousZeroPoleGain = ZeroPoleGain<T, eDiscretization::CONTINUOUS>;
template <class T>
using DiscreteZeroPoleGain = ZeroPoleGain<T, eDiscretization::DISCRETE>;


} // namespace dspbb