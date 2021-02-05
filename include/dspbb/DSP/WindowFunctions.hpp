#pragma once

#include "../Utility/TypeTraits.hpp"
#include "../Primitives/Signal.hpp"

#include <cmath>


namespace enl {

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
Signal<T, Domain> HammingWindow(size_t length) {
	using UnderlyingT = remove_complex_t<T>;

	Signal<T, Domain> window;
	window.Resize(length);
	for (size_t i = 0; i < length; ++i) {
		const UnderlyingT value = UnderlyingT(0.54) - UnderlyingT(0.46) * std::cos((UnderlyingT(2 * 4) * std::atan(UnderlyingT(1)) * i) / UnderlyingT(length - 1));
		window[i] = value;
	}
	return window;
}


template <class T, eSignalDomain Domain = eSignalDomain::TIME>
Signal<T, Domain> KaiserWindow(size_t length, T alpha = T(0.5)) {
	using UnderlyingT = remove_complex_t<T>;

	Signal<T, Domain> window;
	window.Resize(length);
	for (size_t i = 0; i < length; ++i) {
		constexpr UnderlyingT pi = UnderlyingT(3.141592653589793238);
		const UnderlyingT x = 2 * UnderlyingT(i) / UnderlyingT(length) - 1;
		const UnderlyingT value = UnderlyingT(std::cyl_bessel_i(UnderlyingT(0), pi * alpha * std::sqrt(1 - x * x)) / std::cyl_bessel_i(UnderlyingT(0), pi * alpha));
		window[i] = value;
	}
	return window;
}



} // namespace enl