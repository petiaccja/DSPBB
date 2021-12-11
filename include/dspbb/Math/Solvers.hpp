#pragma once

namespace dspbb {

template <class Func, class T>
T Bisect(Func func, T a, T b, T tolerance = T(0)) {
	T c = (a + b) / T(2);
	T err = std::numeric_limits<T>::max();
	T errPrev = err;
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
		errPrev = err;
		err = std::abs(b - a);
	} while (err > tolerance && err < errPrev);

	return c;
}

template <class Func, class Derivative, class T>
T NewtonRaphson(Func func, Derivative der, T x0, T tolerance = T(0)) {
	T x1 = x0;
	T err = std::numeric_limits<T>::max();
	T errPrev = err;
	do {
		x0 = x1;
		x1 = x0 - func(x0) / der(x0);
		errPrev = err;
		err = std::abs(x1 - x0);
	} while (err > tolerance && err < errPrev);
	return x1;
}

} // namespace dspbb