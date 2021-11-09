#pragma once

#include "../../Generators/Spaces.hpp"
#include "../../Utility/Numbers.hpp"
#include "../../LTISystems/Systems.hpp"


namespace dspbb {

template <class T>
ZeroPoleGain<T, eDiscretization::CONTINUOUS> Butterworth(size_t order) {
	FactoredPolynomial<T> poles;
	poles.Resize(order % 2, order / 2, -T(1));

	size_t index = 0;
	for (auto& root : poles.ComplexPairs()) {
		const T phase = pi_v<T> * (T(0.5) + T(index) / T(order) + T(0.5) / T(order));
		root = std::polar(T(1), std::real(phase));
		++index;
	}

	return { T(1), {}, poles };
}


} // namespace dspbb