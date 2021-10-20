#pragma once

#include "../../Generators/Spaces.hpp"
#include "../../Utility/Numbers.hpp"
#include "../../LTISystems/System.hpp"


namespace dspbb {

template <class T>
PoleZeroSystem<T, eSystemDiscretization::CONTINUOUS> Butterworth(unsigned order) {
	FactoredPolynomial<T> poles;
	poles.Resize(order % 2, order / 2, -T(1));
	
	LinSpace(poles.ComplexRoots(), pi_v<T> * T(1.0), pi_v<T> * T(1.5), false);
	poles.ComplexRoots() += pi_v<T> / T(2) / T(order);
	std::for_each(poles.ComplexRoots().begin(), poles.ComplexRoots().end(), [](auto& arg) {
		arg = std::polar(T(1), std::real(arg));
	});

	return { T(1), {}, poles };
}


} // namespace dspbb