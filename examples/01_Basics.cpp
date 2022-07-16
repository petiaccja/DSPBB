//------------------------------------------------------------------------------
// 01. Basics
//
// This example introduces the library and explains the most fundamental classes
// and design philosophy.
//------------------------------------------------------------------------------


// Includes are grouped in different folders. Let's include Signal, the class you
// will see the most, and some other functionality.
#include <dspbb/Filtering/Windowing.hpp>
#include <dspbb/Generators/Waveforms.hpp>
#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Primitives/SignalView.hpp>

#include <iostream>


// Everything that you should use is in the namespace dspbb. There are other,
// potentially useful things in nested namespaces, but they are not meant
// to be stable.
using namespace dspbb;


int main() {
	// When it comes to signal processing, you will mostly work with signals,
	// spectra, or cepstra.
	Signal<float> signal(1024);
	Spectrum<std::complex<float>> spectrum(1024);

	// DSPBB treats signals and spectra as different types so that you don't
	// accidentally mix them up in your code. Arithmetic operators and functions
	// only work with matching types. They are, however, powered by
	// the same underlying implementation and thus have identical interfaces.
	// Here is how they are defined:
	using PlainSignal = BasicSignal<float, DOMAINLESS>;

	// Most DSPBB functions come in two flavors:
	// 1) set an existing memory region, "signal", to contain a square wave,
	SquareWave(signal, 1024, 10.0);
	// 2) return the requested signal in a brand-new memory region.
	const auto window = BlackmanHarrisWindow<float, TIME_DOMAIN>(signal.size());
	// Use the first method when you want to avoid allocation for safety or performance
	// reasons. Otherwise, the second one is often cleaner due to immutability.

	// You can use operators naturally. If you don't want allocations,
	// there are 3-operand function (e.g. Multiply).
	const auto windowed = signal * window;

	// I think you now know where this is going. Note that we could have also
	// used the allocating flavor of the Fft function, but we already had "spectrum" allocated.
	Fft(spectrum, windowed);

	// Signal views can help work with part of signals. In this case, they can be used
	// to verify if the FFT of a real signal indeed contains complex conjugates.
	const SpectrumView<const std::complex<float>> positiveHalf{ spectrum.begin() + 1, 511 };
	const Spectrum<std::complex<float>> negativeHalf{ spectrum.rbegin(), spectrum.rbegin() + 511 };
	const float error = Max(Imag(positiveHalf + negativeHalf));
	std::cout << "Error of FFT's conjugate symmetry: " << error << std::endl;
}