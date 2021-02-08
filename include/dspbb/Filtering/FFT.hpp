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
/// <returns> The complex frequency domain signal, same size as <paramref name="signal"/>. </returns>
template <class T, class Complex = std::enable_if_t<!is_complex_v<T>, std::complex<T>>>
Spectrum<Complex> FourierTransform(TimeSignalView<const T> signal) {
	Spectrum<Complex> fft(signal.Size());
	pocketfft::shape_t shape = { signal.Size() };
	pocketfft::stride_t stride_in = { sizeof(T) };
	pocketfft::stride_t stride_out = { sizeof(Complex) };
	pocketfft::r2c(shape, stride_in, stride_out, 0, pocketfft::FORWARD, signal.Data(), fft.Data(), T(std::sqrt(2.0)));
	return fft;
}

/// <summary>
/// You'd never guess: it computes the discrete fourier transform of a complex <paramref name="signal"/>!
/// </summary>
/// <param name="signal"> Length of the signal is unrestricted. </param>
/// <returns> The complex frequency domain signal, same size as <paramref name="signal"/>. </returns>
template <class T, class Complex = std::enable_if_t<is_complex_v<T>, std::complex<T>>>
Spectrum<T> FourierTransform(TimeSignalView<const T> signal) {
	Spectrum<T> fft(signal.Size());
	pocketfft::shape_t shape = { signal.Size() };
	pocketfft::stride_t stride = { sizeof(T) };
	pocketfft::shape_t axes = { 0 };
	pocketfft::c2c(shape, stride, stride, axes, pocketfft::FORWARD, signal.Data(), fft.Data(), typename T::value_type(std::sqrt(2.0)));
	return fft;
}

/// <summary>
/// Computes the inverse fourier transform of a complex spectrum.
/// </summary>
/// <typeparam name="T"> Please specify T explicitly to an std::complex of your choice. </typeparam>
/// <param name="fft"> Length of the spectrum is unrestricted. </param>
/// <returns> The complex time domain signal, same size as <paramref name="fft"/>. </returns>
template <class T, class Complex = std::enable_if_t<is_complex_v<T>, std::complex<T>>>
TimeSignal<T> FourierTransform(SpectrumView<const T> fft) {
	using U = typename T::value_type;
	TimeSignal<T> signal(fft.Size());
	pocketfft::shape_t shape = { fft.Size() };
	pocketfft::stride_t stride = { sizeof(T) };
	pocketfft::shape_t axes = { 0 };
	pocketfft::c2c(shape, stride, stride, axes, pocketfft::BACKWARD, fft.Data(), signal.Data(), U(1.0 / std::sqrt(2.0) / double(fft.Size())));
	return signal;
}

/// <summary>
/// Computes the inverse fourier transform of a complex spectrum.
/// </summary>
/// <param name="fft"> It better be symmetric. Length of the spectrum is unrestricted. </param>
/// <typeparam name="T"> Please specify T explicitly to a floating type of your choice. </typeparam>
/// <returns> The real time domain signal, same size as <paramref name="fft"/>. </returns>
template <class T, class Complex = std::enable_if_t<!is_complex_v<T>, std::complex<T>>>
TimeSignal<T> FourierTransform(SpectrumView<const std::complex<T>> fft) {
	TimeSignal<T> signal(fft.Size());
	pocketfft::shape_t shape = { fft.Size() };
	pocketfft::stride_t stride_in = { sizeof(std::complex<T>) };
	pocketfft::stride_t stride_out = { sizeof(T) };
	pocketfft::c2r<T>(shape, stride_in, stride_out, 0, pocketfft::BACKWARD, fft.Data(), signal.Data(), T(1.0 / std::sqrt(2.0) / double(fft.Size())));
	return signal;
}


template <class T>
auto FourierTransform(const TimeSignal<T>& signal) {
	return FourierTransform<T>(AsConstView(signal));
}


template <class T, class Complex = std::enable_if_t<is_complex_v<T>, std::complex<T>>>
auto FourierTransform(const Spectrum<T>& fft) {
	return FourierTransform<T>(AsConstView(fft));
}

template <class T, class Complex = std::enable_if_t<!is_complex_v<T>, std::complex<T>>>
auto FourierTransform(const Spectrum<std::complex<T>>& fft) {
	return FourierTransform<T>(AsConstView(fft));
}


inline double FourierBin2Frequency(size_t binIdx, size_t numBins, uint64_t sampleRate) {
	return double(binIdx) / double(numBins) * double(sampleRate);
}

inline size_t FourierFrequency2Bin(double frequency, size_t numBins, uint64_t sampleRate) {
	return size_t(std::round(frequency / double(sampleRate) * double(numBins)));
}


} // namespace dspbb
