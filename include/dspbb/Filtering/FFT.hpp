#pragma once

#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "../Utility/TypeTraits.hpp"

#include <pocketfft/pocketfft_hdronly.h>


namespace dspbb {


/// <summary>
/// You'd never guess: it computes the discrete fourier transform of a real <paramref name="signal"/>!
/// </summary>
/// <param name="signal"> Length of the signal is unrestricted. </param>
/// <param name="full"> If true, the full FFT is returned, its length the same as the signal.
///		If false, only the first (almost) half is returned, the rest can be derived from that as
///		the DFT is even-symmetric for real inputs. </param>
/// <returns> The complex frequency domain signal, same size as <paramref name="signal"/>. </returns>
template <class T, class Complex = std::enable_if_t<!is_complex_v<T>, std::complex<T>>>
Spectrum<Complex> FourierTransform(TimeSignalView<const T> signal, bool full = false) {
	const size_t outputSize = full ? signal.Size() : signal.Size() / 2 + 1;
	Spectrum<Complex> fft(outputSize);
	pocketfft::shape_t shape = { signal.Size() };
	pocketfft::stride_t stride_in = { sizeof(T) };
	pocketfft::stride_t stride_out = { sizeof(Complex) };
	pocketfft::r2c(shape, stride_in, stride_out, 0, pocketfft::FORWARD, signal.Data(), fft.Data(), T(1));
	if (full) {
		for (size_t i = 1; i < (outputSize + 1) / 2; ++i) {
			fft[outputSize - i] = std::conj(fft[i]);
		}
	}
	return fft;
}

/// <summary>
/// You'd never guess: it computes the discrete fourier transform of a complex <paramref name="signal"/>!
/// </summary>
/// <param name="signal"> Length of the signal is unrestricted. </param>
/// <returns> The complex frequency domain signal, same size as <paramref name="signal"/>. </returns>
template <class T>
Spectrum<std::complex<T>> FourierTransform(TimeSignalView<const std::complex<T>> signal) {
	Spectrum<std::complex<T>> fft(signal.Size());
	pocketfft::shape_t shape = { signal.Size() };
	pocketfft::stride_t stride = { sizeof(std::complex<T>) };
	pocketfft::shape_t axes = { 0 };
	pocketfft::c2c(shape, stride, stride, axes, pocketfft::FORWARD, signal.Data(), fft.Data(), T(1));
	return fft;
}

/// <summary>
/// Computes the inverse fourier transform of a complex spectrum.
/// </summary>
/// <typeparam name="T"> Please specify T explicitly to an std::complex of your choice. </typeparam>
/// <param name="fft"> Length of the spectrum is unrestricted. </param>
/// <returns> The complex time domain signal, same size as <paramref name="fft"/>. </returns>
template <class T>
TimeSignal<std::complex<T>> InverseFourierTransformC(SpectrumView<const std::complex<T>> fft) {
	TimeSignal<std::complex<T>> signal(fft.Size());
	pocketfft::shape_t shape = { fft.Size() };
	pocketfft::stride_t stride = { sizeof(std::complex<T>) };
	pocketfft::shape_t axes = { 0 };
	pocketfft::c2c(shape, stride, stride, axes, pocketfft::BACKWARD, fft.Data(), signal.Data(), T(1.0 / double(fft.Size())));
	return signal;
}

/// <summary>
/// Computes the inverse fourier transform of a complex spectrum.
/// </summary>
/// <param name="fft"> If <paramref name="full"/>, only the first half is taken, and the output is
///		the same length as the <paramref name="fft"/>. If not <paramref name="full"/>,
///		then the whole is taken, and the output is (about) twice the length of the <paramref name="fft"/>.
///		</param>
/// <param name="full"> Tells if <paramref name="fft"/> contains the whole, symmetric FFT,
///		or only the non-redundant first half. </param>
/// <typeparam name="T"> Please specify T explicitly to a floating type of your choice. </typeparam>
/// <returns> The real time domain signal. </returns>
template <class T>
TimeSignal<T> InverseFourierTransformR(SpectrumView<const std::complex<T>> fft, bool full = false) {
	assert(!fft.Empty());
	const size_t signalSize = full ? fft.Size() : fft.Size() * 2 - 1 - fft.Size() % 2;
	TimeSignal<T> signal(signalSize);
	pocketfft::shape_t shape = { signalSize };
	pocketfft::stride_t stride_in = { sizeof(std::complex<T>) };
	pocketfft::stride_t stride_out = { sizeof(T) };
	pocketfft::c2r<T>(shape, stride_in, stride_out, 0, pocketfft::BACKWARD, fft.Data(), signal.Data(), T(1.0 / double(signalSize)));
	return signal;
}



template <class T, class Complex = std::enable_if_t<!is_complex_v<T>, std::complex<T>>>
Spectrum<Complex> FourierTransform(const TimeSignal<T>& signal, bool full = false) {
	return FourierTransform(AsConstView(signal), full);
}

template <class T>
auto FourierTransform(const TimeSignal<std::complex<T>>& signal) {
	return FourierTransform(AsConstView(signal));
}

template <class T>
auto InverseFourierTransformC(const Spectrum<std::complex<T>>& fft) {
	return InverseFourierTransformC(AsConstView(fft));
}

template <class T>
auto InverseFourierTransformR(const Spectrum<std::complex<T>>& fft, bool full = false) {
	return InverseFourierTransformR(AsConstView(fft), full);
}


inline double FourierBin2Frequency(size_t binIdx, size_t numBins, uint64_t sampleRate) {
	return double(binIdx) / double(numBins) * double(sampleRate);
}

inline size_t FourierFrequency2Bin(double frequency, size_t numBins, uint64_t sampleRate) {
	return size_t(std::round(frequency / double(sampleRate) * double(numBins)));
}


} // namespace dspbb
