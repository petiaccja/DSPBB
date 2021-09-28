#pragma once

#include "../../Generators/Spaces.hpp"
#include "../../Utility/Numbers.hpp"
#include "PoleZeroSystem.hpp"


namespace dspbb {

template <class T>
PoleZeroSystem<T, eSystemDiscretization::CONTINUOUS> Butterworth(unsigned order) {
	std::vector<std::complex<T>> poles(order);

	// Calculate evenly spaced poles on real+ imaginary- quarter plane.
	const unsigned numPositiveImaginaryPoles = order / 2;
	SignalView<std::complex<T>, eSignalDomain::DOMAINLESS> polesView{ poles.data(), order };
	LinSpace(polesView, pi_v<T> * T(0.5), pi_v<T> * T(1.5), false);
	polesView += pi_v<T> / T(2) / T(order);
	std::for_each(poles.begin(), poles.begin() + numPositiveImaginaryPoles, [](auto& arg) {
		arg = std::polar(T(1), std::real(arg));
	});

	// Mirror-copy poles across real axis.
	std::transform(poles.begin(), poles.begin() + numPositiveImaginaryPoles, poles.rbegin(), [](const auto& arg) {
		return std::conj(arg);
	});

	// Add lone pole on the real axis.
	if (order % 2 == 1) {
		poles[order / 2] = -T(1);
	}

	return { T(1), poles, {} };
}


} // namespace dspbb