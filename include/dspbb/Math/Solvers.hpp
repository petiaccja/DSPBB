#pragma once

namespace dspbb {

template <class Func, class T>
T Bisect(Func func, T a, T b, T tolerance = T(0)) {
	if (b < a) {
		std::swap(a, b);
	}
	T c = (a + b) / T(2);
	do {
		const auto fa = func(a);
		const auto fc = func(c);
		if (fa * fc <= 0) {
			b = c;
		}
		else {
			a = c;
		}
		c = (a + b) / T(2);
	} while (b - c > tolerance && c - a > tolerance);

	return c;
}

template <class Func, class Derivative, class T>
T NewtonRaphson(Func func, Derivative der, T x0, T tolerance = T(0)) {
	T x1 = x0;
	do {
		x0 = x1;
		x1 = x0 - func(x0) / der(x0);
	} while (std::abs(x1 - x0) > tolerance);
	return x1;
}

} // namespace dspbb