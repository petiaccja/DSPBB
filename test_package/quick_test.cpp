//L=============================================================================
//L This software is distributed under the MIT license.
//L Copyright 2021 Péter Kardos
//L=============================================================================

#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Filtering/FFT.hpp>
#include <iostream>

int main() {
    const dspbb::TimeSignalF u = { 1, 2, 3 };
    const dspbb::TimeSignalF v = { 15, 23, 32 };
	const auto c = u + v;
	const dspbb::SpectrumCF s = FourierTransform(u);
	std::cout << "DSPBB installation works." << std::endl;
	return 0;
}
