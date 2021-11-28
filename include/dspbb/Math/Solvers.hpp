#pragma once

namespace dspbb {

template <class Func, class T>
T Bisect(Func func, T a, T b) {
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
	} while (c != a && c != b);

	return c;
}

} // namespace dspbb