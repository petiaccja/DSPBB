#pragma once

#include "dspbb/Primitives/Signal.hpp"
#include "dspbb/Primitives/SignalView.hpp"

#include <type_traits>


namespace dspbb {


template <class T, class U>
using ProductT = decltype(std::declval<T>() * std::declval<U>());


template <class T, class U, eSignalDomain Domain>
ProductT<T, U> DotProduct(SignalView<const T, Domain> s1, SignalView<const U, Domain> s2, size_t length) {
	ProductT<T, U> sum{ 0 };
	for (size_t i = 0; i < length; ++i) {
		sum += s1[i] * s2[i];
	}
	return sum;
}

template <class T, class U, eSignalDomain Domain>
std::complex<ProductT<T, U>> DotProduct(SignalView<const std::complex<T>, Domain> s1, SignalView<const U, Domain> s2, size_t length) {
	using OutputT = std::complex<ProductT<T, U>>;
	OutputT sum{ 0 };
	for (size_t i = 0; i < length; ++i) {
		sum += OutputT(s1[i]) * OutputT(s2[i]);
	}
	return sum;
}

template <class T, class U, eSignalDomain Domain>
std::complex<ProductT<T, U>> DotProduct(SignalView<const T, Domain> s1, SignalView<const std::complex<U>, Domain>& s2, size_t length) {
	using OutputT = std::complex<ProductT<T, U>>;
	OutputT sum{ 0 };
	for (size_t i = 0; i < length; ++i) {
		sum += OutputT(s1[i]) * OutputT(s2[i]);
	}
	return sum;
}

template <class T, class U, eSignalDomain Domain>
std::complex<ProductT<T, U>> DotProduct(SignalView<const std::complex<T>, Domain> s1, SignalView<const std::complex<U>, Domain> s2, size_t length) {
	using OutputT = std::complex<ProductT<T, U>>;
	OutputT sum{ 0 };
	for (size_t i = 0; i < length; ++i) {
		sum += OutputT(s1[i]) * OutputT(s2[i]);
	}
	return sum;
}

template <class T, class U, eSignalDomain Domain>
auto DotProduct(const Signal<T, Domain>& s1, SignalView<const U, Domain> s2, size_t length) {
	return DotProduct(AsConstSpan(s1), s2, length);
}

template <class T, class U, eSignalDomain Domain>
auto DotProduct(SignalView<const T, Domain> s1, const Signal<U, Domain>& s2, size_t length) {
	return DotProduct(s1, AsConstSpan(s2), length);
}

template <class T, class U, eSignalDomain Domain>
auto DotProduct(const Signal<T, Domain>& s1, const Signal<U, Domain>& s2, size_t length) {
	return DotProduct(AsConstSpan(s1), AsConstSpan(s2), length);
}


} // namespace dspbb