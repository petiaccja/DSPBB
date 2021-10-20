#pragma once

#include "../../Generators/Spaces.hpp"
#include "../../Utility/Numbers.hpp"
#include "../../LTISystems/System.hpp"


namespace dspbb {

template <class T>
PoleZeroSystem<T, eSystemDiscretization::CONTINUOUS> Butterworth(unsigned order) {
	FactoredPolynomial<T> poles;
	poles.Resize(order % 2, order / 2, -T(1));

	size_t index = 0;
	for (auto& root : poles.ComplexRoots()) {
		const T phase = pi_v<T> * (T(0.5) + T(index) / T(order) + T(0.5) / T(order));
		root = std::polar(T(1), std::real(phase));
		++index;
	}

	return { T(1), {}, poles };
}


} // namespace dspbb