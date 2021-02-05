#pragma once

#include "../Primitives/Signal.hpp"
#include "../Utility/TypeTraits.hpp"

#include <pocketfft/pocketfft_hdronly.h>


namespace dspbb {


template <class T, class Complex = std::enable_if_t<!is_complex_v<T>, std::complex<T>>>
Spectrum<Complex> FourierTransform(const TimeSignal<T>& signal) {
	Spectrum<Complex> fft(signal.Size());
	pocketfft::shape_t shape = { signal.Size() };
	pocketfft::stride_t stride_in = { sizeof(T) };
	pocketfft::stride_t stride_out = { sizeof(Complex) };
	pocketfft::r2c(shape, stride_in, stride_out, 0, pocketfft::FORWARD, signal.Data(), fft.Data(), T(std::sqrt(2.0)));
	return fft;
}

template <class T, class Complex = std::enable_if_t<is_complex_v<T>, std::complex<T>>>
Spectrum<T> FourierTransform(const TimeSignal<T>& signal) {
	Spectrum<T> fft(signal.Size());
	pocketfft::shape_t shape = { signal.Size() };
	pocketfft::stride_t stride = { sizeof(T) };
	pocketfft::shape_t axes = { 0 };
	pocketfft::c2c(shape, stride, stride, axes, pocketfft::FORWARD, signal.Data(), fft.Data(), T::value_type(std::sqrt(2.0)));
	return fft;
}


template <class T, class Complex = std::enable_if_t<is_complex_v<T>, std::complex<T>>>
TimeSignal<T> FourierTransform(const Spectrum<T>& fft) {
	using U = typename T::value_type;
	TimeSignal<T> signal(fft.Size());
	pocketfft::shape_t shape = { fft.Size() };
	pocketfft::stride_t stride = { sizeof(T) };
	pocketfft::shape_t axes = { 0 };
	pocketfft::c2c(shape, stride, stride, axes, pocketfft::BACKWARD, fft.Data(), signal.Data(), U(1.0 / std::sqrt(2.0) / double(fft.Size())));
	return signal;
}

template <class T, class Complex = std::enable_if_t<!is_complex_v<T>, std::complex<T>>>
TimeSignal<T> FourierTransform(const Spectrum<std::complex<T>>& fft) {
	TimeSignal<T> signal(fft.Size());
	pocketfft::shape_t shape = { fft.Size() };
	pocketfft::stride_t stride_in = { sizeof(std::complex<T>) };
	pocketfft::stride_t stride_out = { sizeof(T) };
	pocketfft::c2r<T>(shape, stride_in, stride_out, 0, pocketfft::BACKWARD, fft.Data(), signal.Data(), T(1.0 / std::sqrt(2.0) / double(fft.Size())));
	return signal;
}

inline double FourierBin2Frequency(size_t binIdx, size_t numBins, uint64_t sampleRate) {
	return double(binIdx) / double(numBins) * double(sampleRate);
}

inline size_t FourierFrequency2Bin(double frequency, size_t numBins, uint64_t sampleRate) {
	return size_t(std::round(frequency / double(sampleRate) * double(numBins)));
}


} // namespace dspbb
