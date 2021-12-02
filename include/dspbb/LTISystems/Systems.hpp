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

	size_t Order() const { return std::max(zeros.Order(), poles.Order()); }
};


template <class T, eDiscretization Discretization>
struct TransferFunction {
	TransferFunction() = default;
	explicit TransferFunction(const ZeroPoleGain<T, Discretization>& zpk);

	Polynomial<T> numerator;
	Polynomial<T> denominator;

	std::complex<T> operator()(const std::complex<T>& x) const { return numerator(x) / denominator(x); }
	T operator()(const T& x) const { return numerator(x) / denominator(x); }

	size_t Order() const { return std::max(numerator.Order(), denominator.Order()); }
};


template <class T>
struct CascadedBiquad {
	CascadedBiquad() = default;
	explicit CascadedBiquad(const ZeroPoleGain<T, eDiscretization::DISCRETE>& zpk);

	struct Biquad {
		std::array<T, 3> numerator = { 0, 0, 1 };
		std::array<T, 2> denominator = { 0, 0 };
		uint8_t numOrder = 0;
		uint8_t denOrder = 0;
	};
	std::vector<Biquad> sections;

	std::complex<T> operator()(const std::complex<T>& x) const;
	T operator()(const T& x) const;

	size_t Order() const;

private:
	template <class Iter>
	static std::vector<std::pair<uint8_t, std::array<T, 3>>> RealRootPolynomials(Iter first, Iter last);
	template <class Iter>
	static std::vector<std::pair<uint8_t, std::array<T, 3>>> ComplexPairPolynomials(Iter first, Iter last);
	template <class X>
	static X EvalSection(const Biquad& section, const X& x);
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
	zeroPolys.insert(zeroPolys.end(), complexZeroPolys.begin(), complexZeroPolys.end());
	auto polePolys = realPolePolys;
	polePolys.insert(polePolys.end(), complexPolePolys.begin(), complexPolePolys.end());

	const auto polyAscending = [](const auto& lhs, const auto& rhs) {
		return std::abs(lhs.second[0]) < std::abs(rhs.second[0]);
	};
	std::sort(zeroPolys.begin(), zeroPolys.end(), polyAscending);
	std::sort(polePolys.begin(), polePolys.end(), polyAscending);

	const size_t numSections = std::max(zeroPolys.size(), polePolys.size());
	sections.resize(numSections);
	for (size_t i = 0; i < zeroPolys.size(); ++i) {
		sections[i].numerator = zeroPolys[i].second;
		sections[i].numOrder = zeroPolys[i].first;
	}
	for (size_t i = 0; i < polePolys.size(); ++i) {
		sections[i].denominator = { polePolys[i].second[0], polePolys[i].second[1] };
		sections[i].denOrder = polePolys[i].first;
	}
	if (!sections.empty()) {
		sections.back().numerator[0] *= zpk.gain;
		sections.back().numerator[1] *= zpk.gain;
		sections.back().numerator[2] *= zpk.gain;
	}
}

template <class T>
std::complex<T> CascadedBiquad<T>::operator()(const std::complex<T>& x) const {
	return std::transform_reduce(sections.begin(), sections.end(), std::complex<T>(T(1)), std::multiplies{}, [&x](const auto& section) {
		return EvalSection(section, x);
	});
}

template <class T>
T CascadedBiquad<T>::operator()(const T& x) const {
	return std::transform_reduce(sections.begin(), sections.end(), T(1), std::multiplies{}, [&x](const auto& section) {
		return EvalSection(section, x);
	});
}

template <class T>
size_t CascadedBiquad<T>::Order() const {
	size_t numOrder = 1;
	size_t denOrder = 1;
	for (auto& section : sections) {
		numOrder += section.numOrder;
		denOrder += section.denOrder;
	}
	return sections.empty() ? 0 : std::max(numOrder, denOrder);
}

template <class T>
template <class Iter>
std::vector<std::pair<uint8_t, std::array<T, 3>>> CascadedBiquad<T>::RealRootPolynomials(Iter first, Iter last) {
	std::vector<T> ascending{ first, last };
	std::sort(ascending.begin(), ascending.end());

	std::vector<std::pair<uint8_t, std::array<T, 3>>> polynomials;

	const auto partitionPoint = std::partition_point(ascending.begin(), ascending.end(), [](const T& r) { return r < T(0); });
	auto firstNegative = ascending.begin();
	auto lastNegative = std::reverse_iterator{ partitionPoint };
	auto firstPositive = partitionPoint;
	auto lastPositive = std::reverse_iterator{ ascending.end() };
	for (; lastNegative.base() - firstNegative > 1; ++firstNegative, ++lastNegative) {
		const auto& r1 = *firstNegative;
		const auto& r2 = *lastNegative;
		polynomials.emplace_back(2, std::array{ r1 * r2, -r1 - r2, T(1) });
	}
	for (; lastPositive.base() - firstPositive > 1; ++firstPositive, ++lastPositive) {
		const auto& r1 = *firstPositive;
		const auto& r2 = *lastPositive;
		polynomials.emplace_back(2, std::array{ r1 * r2, -r1 - r2, T(1) });
	}
	const bool unpairedNegative = firstNegative < lastNegative.base();
	const bool unpairedPositive = firstPositive < lastPositive.base();
	if (unpairedNegative && unpairedPositive) {
		const auto& r1 = *firstPositive;
		const auto& r2 = *firstNegative;
		polynomials.emplace_back(2, std::array{ r1 * r2, -r1 - r2, T(1) });
	}
	else if (unpairedPositive) {
		polynomials.emplace_back(1, std::array{ T(0), -*firstPositive, T(1) });
	}
	else if (unpairedNegative) {
		polynomials.emplace_back(1, std::array{ T(0), -*firstNegative, T(1) });
	}

	return polynomials;
}

template <class T>
template <class Iter>
std::vector<std::pair<uint8_t, std::array<T, 3>>> CascadedBiquad<T>::ComplexPairPolynomials(Iter first, Iter last) {
	std::vector<std::pair<uint8_t, std::array<T, 3>>> polynomials;
	while (first != last) {
		std::complex<T> r = *first;
		polynomials.emplace_back(2, std::array{ std::abs(r) * std::abs(r), -T(2) * std::real(r), T(1) });
		++first;
	}
	return polynomials;
}

template <class T>
template <class X>
X CascadedBiquad<T>::EvalSection(const Biquad& section, const X& x) {
	const std::array xs = { X(T(0)),
							X(T(1)),
							x,
							x * x };
	const auto num = section.numerator[0] + xs[section.numOrder] * section.numerator[1] + xs[1 + section.numOrder] * section.numerator[2];
	const auto den = section.denominator[0] + xs[section.denOrder] * section.denominator[1] + xs[1 + section.denOrder];
	return num / den;
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