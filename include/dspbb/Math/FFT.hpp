#pragma once

#include "../PocketFFT/pocketfft_hdronly.h"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "../Utility/TypeTraits.hpp"

#include <algorithm>


namespace dspbb {

namespace impl {
	template <class T>
	void Fft(SpectrumView<std::complex<T>> out, TimeSignalView<const T> in) {
		const size_t halfSize = in.Size() / 2 + 1;
		const size_t fullSize = in.Size();
		if (out.Size() != halfSize || out.Size() != fullSize) {
			throw std::invalid_argument("Output size must be N or N/2+1.");
		}

		pocketfft_dspbb::shape_t shape = { in.Size() };
		pocketfft_dspbb::stride_t stride_in = { sizeof(T) };
		pocketfft_dspbb::stride_t stride_out = { sizeof(std::complex<T>) };
		pocketfft_dspbb::r2c(shape, stride_in, stride_out, 0, pocketfft_dspbb::FORWARD, in.Data(), out.Data(), T(1));

		if (out.Size() == fullSize && fullSize > 2) {
			std::reverse_copy(out.begin() + 1, out.begin() + fullSize / 2, out.begin() + fullSize / 2 + 1);
		}
	}

	template <class T>
	void Fft(SpectrumView<std::complex<T>> out, TimeSignalView<const std::complex<T>> in) {
		if (out.Size() != in.Size()) {
			throw std::invalid_argument("Output and input size must be the same.");
		}

		pocketfft_dspbb::shape_t shape = { in.Size() };
		pocketfft_dspbb::stride_t stride = { sizeof(std::complex<T>) };
		pocketfft_dspbb::shape_t axes = { 0 };
		pocketfft_dspbb::c2c(shape, stride, stride, axes, pocketfft_dspbb::FORWARD, in.Data(), out.Data(), T(1));
	}

	template <class T>
	void Ifft(TimeSignalView<T> out, TimeSignalView<const std::complex<T>> in) {
		const size_t halfSize = out.Size() / 2 + 1;
		const size_t fullSize = out.Size();
		if (in.Size() != halfSize || in.Size() != fullSize) {
			throw std::invalid_argument("Input size must be N or N/2+1.");
		}

		pocketfft_dspbb::shape_t shape = { out.Size() };
		pocketfft_dspbb::stride_t stride_in = { sizeof(std::complex<T>) };
		pocketfft_dspbb::stride_t stride_out = { sizeof(T) };
		pocketfft_dspbb::c2r<T>(shape, stride_in, stride_out, 0, pocketfft_dspbb::BACKWARD, in.Data(), out.Data(), T(1.0 / double(out.Size())));
	}

	template <class T>
	void Ifft(TimeSignalView<std::complex<T>> out, TimeSignalView<const std::complex<T>> in) {
		if (out.Size() != in.Size()) {
			throw std::invalid_argument("Output and input size must be the same.");
		}

		pocketfft_dspbb::shape_t shape = { out.Size() };
		pocketfft_dspbb::stride_t stride = { sizeof(std::complex<T>) };
		pocketfft_dspbb::shape_t axes = { 0 };
		pocketfft_dspbb::c2c(shape, stride, stride, axes, pocketfft_dspbb::BACKWARD, in.Data(), out.Data(), T(1.0 / double(out.Size())));
	}


	template <class T>
	Spectrum<std::complex<T>> Fft(TimeSignalView<const T> in, bool full) {
		const size_t halfSize = in.Size() / 2 + 1;
		const size_t fullSize = in.Size();
		const size_t size = full ? fullSize : halfSize;

		Spectrum<std::complex<T>> out(size);
		Fft(out, in);
		return out;
	}

	template <class T>
	Spectrum<std::complex<T>> Fft(TimeSignalView<const std::complex<T>> in) {
		const size_t size = in.Size();

		Spectrum<std::complex<T>> out(size);
		Fft(out, in);
		return out;
	}

	template <class T>
	TimeSignal<T> Ifft(SpectrumView<const std::complex<T>> in, bool full, bool even) {
		const size_t halfSizeEven = in.Size() * 2 - 2;
		const size_t halfSizeOdd = in.Size() * 2 - 1;
		const size_t fullSize = in.Size();
		const size_t size = full ? fullSize : (even ? halfSizeEven : halfSizeOdd);

		TimeSignal<std::complex<T>> out(size);
		Fft(out, in);
		return out;
	}

	template <class T>
	TimeSignal<std::complex<T>> Ifft(SpectrumView<const std::complex<T>> in) {
		const size_t size = in.Size();

		TimeSignal<std::complex<T>> out(size);
		Fft(out, in);
		return out;
	}

} // namespace impl


template <class SignalR, class SignalT>
auto Fft(SignalR&& out, const SignalT& in) -> decltype(impl::Fft(AsView(out), AsView(in))) {
	return impl::Fft(AsView(out), AsView(in));
}

template <class SignalR, class SignalT>
auto Ifft(SignalR&& out, const SignalT& in) -> decltype(impl::Fft(AsView(out), AsView(in))) {
	return impl::Fft(AsView(out), AsView(in));
}


template <class SignalT>
auto Fft(const SignalT& in, bool full) -> decltype(impl::Fft(AsView(in), full)) {
	return impl::Fft(AsView(in, full));
}

template <class SignalT>
auto Fft(const SignalT& in) -> decltype(impl::Fft(AsView(in))) {
	return impl::Fft(AsView(in));
}

template <class SignalT>
auto Ifft(const SignalT& in, bool full, bool even) -> decltype(impl::Ifft(AsView(in), full, even)) {
	return impl::Ifft(AsView(in), full, even);
}

template <class SignalT>
auto Ifft(const SignalT& in) -> decltype(impl::Fft(AsView(in))) {
	return impl::Fft(AsView(in));
}

/// <summary>
/// You'd never guess: it computes the discrete fourier transform of a real <paramref name="signal"/>!
/// </summary>
/// <param name="signal"> Length of the signal is unrestricted. </param>
/// <param name="full"> If true, the full FFT is returned, its length the same as the signal.
///		If false, only the first (almost) half is returned, the rest can be derived from that as
///		the DFT is even-symmetric for real inputs. </param>
/// <returns> The complex frequency domain signal, same size as <paramref name="signal"/>. </returns>
template <class T, class Complex = std::enable_if_t<!is_complex_v<T>, std::complex<T>>>
Spectrum<Complex> FourierTransform(TimeSignalView<const T> signal, bool full = true) {
	const size_t outputSize = full ? signal.Size() : signal.Size() / 2 + 1;
	Spectrum<Complex> fft(outputSize);
	pocketfft_dspbb::shape_t shape = { signal.Size() };
	pocketfft_dspbb::stride_t stride_in = { sizeof(T) };
	pocketfft_dspbb::stride_t stride_out = { sizeof(Complex) };
	pocketfft_dspbb::r2c(shape, stride_in, stride_out, 0, pocketfft_dspbb::FORWARD, signal.Data(), fft.Data(), T(1));
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
	pocketfft_dspbb::shape_t shape = { signal.Size() };
	pocketfft_dspbb::stride_t stride = { sizeof(std::complex<T>) };
	pocketfft_dspbb::shape_t axes = { 0 };
	pocketfft_dspbb::c2c(shape, stride, stride, axes, pocketfft_dspbb::FORWARD, signal.Data(), fft.Data(), T(1));
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
	pocketfft_dspbb::shape_t shape = { fft.Size() };
	pocketfft_dspbb::stride_t stride = { sizeof(std::complex<T>) };
	pocketfft_dspbb::shape_t axes = { 0 };
	pocketfft_dspbb::c2c(shape, stride, stride, axes, pocketfft_dspbb::BACKWARD, fft.Data(), signal.Data(), T(1.0 / double(fft.Size())));
	return signal;
}

/// <summary>
/// Computes the inverse fourier transform of a complex spectrum.
/// </summary>
/// <param name="fft"> If <paramref name="full"/>, only the first half is taken, and the output is
///		the same length as the <paramref name="fft"/>. If not <paramref name="full"/>,
///		then the whole is taken, and the output is (about) twice the length of the <paramref name="fft"/>.
///		</param>
/// <param name="size"> Tells if <paramref name="fft"/> contains the whole, symmetric FFT,
///		or only the non-redundant first half. If size is zero, the FFT is assumed to be full,
///		otherwise size specifies the size of the original signal. </param>
/// <typeparam name="T"> Please specify T explicitly to a floating type of your choice. </typeparam>
/// <returns> The real time domain signal. </returns>
template <class T>
TimeSignal<T> InverseFourierTransformR(SpectrumView<const std::complex<T>> fft, size_t size = 0) {
	assert(!fft.Empty());
	const size_t signalSize = size == 0 ? fft.Size() : size;
	TimeSignal<T> signal(signalSize);
	pocketfft_dspbb::shape_t shape = { signalSize };
	pocketfft_dspbb::stride_t stride_in = { sizeof(std::complex<T>) };
	pocketfft_dspbb::stride_t stride_out = { sizeof(T) };
	pocketfft_dspbb::c2r<T>(shape, stride_in, stride_out, 0, pocketfft_dspbb::BACKWARD, fft.Data(), signal.Data(), T(1.0 / double(signalSize)));
	return signal;
}



template <class T, class Complex = std::enable_if_t<!is_complex_v<T>, std::complex<T>>>
Spectrum<Complex> FourierTransform(const TimeSignal<T>& signal, bool full = true) {
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
auto InverseFourierTransformR(const Spectrum<std::complex<T>>& fft, size_t size = 0) {
	return InverseFourierTransformR(AsConstView(fft), size);
}

template <class T, class Integral, std::enable_if_t<std::is_integral_v<Integral>, int> = 0>
TimeSignal<T> InverseFourierTransformR(SpectrumView<const std::complex<T>> fft, Integral size) {
	static_assert(!std::is_same_v<std::decay_t<Integral>, bool>, "Don't use bool here, that overload has been replaced.");
	return InverseFourierTransformR(fft, (size_t)size);
}

template <class T, class Integral, std::enable_if_t<std::is_integral_v<Integral>, int> = 0>
auto InverseFourierTransformR(const Spectrum<std::complex<T>>& fft, Integral size) {
	static_assert(!std::is_same_v<std::decay_t<Integral>, bool>, "Don't use bool here, that overload has been replaced.");
	return InverseFourierTransformR(fft, (size_t)size);
}


inline double FourierBin2Frequency(size_t binIdx, size_t numBins, uint64_t sampleRate) {
	return double(binIdx) / double(numBins) * double(sampleRate);
}

inline size_t FourierFrequency2Bin(double frequency, size_t numBins, uint64_t sampleRate) {
	return size_t(std::round(frequency / double(sampleRate) * double(numBins)));
}

template <class SignalR, class SignalT, std::enable_if_t<is_same_domain_v<SignalR, SignalT>, int> = 0>
void FftShift(SignalR&& out, const SignalT& in) {
	assert(out.Size() == in.Size());
	if (static_cast<const void*>(std::addressof(out)) == static_cast<const void*>(std::addressof(in))) {
		const auto first = out.begin();
		const auto mid = out.begin() + (1 + out.Size()) / 2;
		const auto last = out.end();
		std::rotate(first, mid, last);
	}
	else {
		const auto first = in.begin();
		const auto mid = in.begin() + (1 + in.Size()) / 2;
		const auto last = in.end();
		const auto write = out.begin();
		std::rotate_copy(first, mid, last, write);
	}
}

template <class SignalT, std::enable_if_t<is_signal_like_v<SignalT>, int> = 0>
SignalT FftShift(const SignalT& in) {
	SignalT out(in.Size());
	FftShift(out, in);
	return out;
}


} // namespace dspbb
