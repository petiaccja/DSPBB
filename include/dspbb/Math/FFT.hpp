#pragma once

#include "../Math/Functions.hpp"
#include "../PocketFFT/pocketfft_hdronly.h"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"

#include <algorithm>


namespace dspbb {


//------------------------------------------------------------------------------
// Kernels
//------------------------------------------------------------------------------
namespace impl {
	struct FftFull {};
	struct FftHalf {};
	constexpr FftFull FFT_FULL;
	constexpr FftHalf FFT_HALF;

	template <class T>
	void Fft(SpectrumView<std::complex<T>> out, TimeSignalView<const T> in) {
		const size_t halfSize = in.Size() / 2 + 1;
		const size_t fullSize = in.Size();
		if (out.Size() != halfSize && out.Size() != fullSize) {
			throw std::invalid_argument("Output size must be N or N/2+1.");
		}

		pocketfft_dspbb::shape_t shape = { in.Size() };
		pocketfft_dspbb::stride_t stride_in = { sizeof(T) };
		pocketfft_dspbb::stride_t stride_out = { sizeof(std::complex<T>) };
		pocketfft_dspbb::r2c(shape, stride_in, stride_out, 0, pocketfft_dspbb::FORWARD, in.Data(), out.Data(), T(1));

		if (out.Size() == fullSize && fullSize > 2) {
			auto first = out.begin() + 1;
			auto last = out.begin() + (fullSize + 1) / 2;
			auto dest = out.begin() + fullSize / 2 + 1;
			std::reverse_copy(first, last, dest);
			const auto mirrorRange = AsView<FREQUENCY_DOMAIN>(dest, out.end());
			Conj(mirrorRange, mirrorRange);
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
	void Ifft(TimeSignalView<T> out, SpectrumView<const std::complex<T>> in) {
		const size_t halfSize = out.Size() / 2 + 1;
		const size_t fullSize = out.Size();
		if (in.Size() != halfSize && in.Size() != fullSize) {
			throw std::invalid_argument("Input size must be N or N/2+1.");
		}

		pocketfft_dspbb::shape_t shape = { out.Size() };
		pocketfft_dspbb::stride_t stride_in = { sizeof(std::complex<T>) };
		pocketfft_dspbb::stride_t stride_out = { sizeof(T) };
		pocketfft_dspbb::c2r<T>(shape, stride_in, stride_out, 0, pocketfft_dspbb::BACKWARD, in.Data(), out.Data(), T(1.0 / double(out.Size())));
	}

	template <class T>
	void Ifft(TimeSignalView<std::complex<T>> out, SpectrumView<const std::complex<T>> in) {
		if (out.Size() != in.Size()) {
			throw std::invalid_argument("Output and input size must be the same.");
		}

		pocketfft_dspbb::shape_t shape = { out.Size() };
		pocketfft_dspbb::stride_t stride = { sizeof(std::complex<T>) };
		pocketfft_dspbb::shape_t axes = { 0 };
		pocketfft_dspbb::c2c(shape, stride, stride, axes, pocketfft_dspbb::BACKWARD, in.Data(), out.Data(), T(1.0 / double(out.Size())));
	}


	template <class T>
	Spectrum<std::complex<T>> Fft(TimeSignalView<const T> in, FftFull) {
		const size_t fullSize = in.Size();
		Spectrum<std::complex<T>> out(fullSize);
		Fft(AsView(out), in);
		return out;
	}

	template <class T>
	Spectrum<std::complex<T>> Fft(TimeSignalView<const T> in, FftHalf) {
		const size_t halfSize = in.Size() / 2 + 1;
		Spectrum<std::complex<T>> out(halfSize);
		Fft(AsView(out), in);
		return out;
	}

	template <class T>
	Spectrum<std::complex<T>> Fft(TimeSignalView<const std::complex<T>> in) {
		const size_t size = in.Size();

		Spectrum<std::complex<T>> out(size);
		Fft(AsView(out), in);
		return out;
	}

	template <class T>
	Signal<T> Ifft(SpectrumView<const std::complex<T>> in, FftHalf, bool even) {
		const size_t halfSizeEven = in.Size() * 2 - 2;
		const size_t halfSizeOdd = in.Size() * 2 - 1;
		Signal<T> out(even ? halfSizeEven : halfSizeOdd);
		Ifft(AsView(out), in);
		return out;
	}

	template <class T>
	Signal<T> Ifft(SpectrumView<const std::complex<T>> in, FftFull) {
		const size_t fullSize = in.Size();
		Signal<T> out(fullSize);
		Ifft(AsView(out), in);
		return out;
	}

	template <class T>
	Signal<std::complex<T>> Ifft(SpectrumView<const std::complex<T>> in) {
		const size_t size = in.Size();
		Signal<std::complex<T>> out(size);
		Ifft(AsView(out), in);
		return out;
	}

} // namespace impl


//------------------------------------------------------------------------------
// Wrappers
//------------------------------------------------------------------------------

using impl::FFT_FULL;
using impl::FFT_HALF;

template <class SignalR, class SignalT>
auto Fft(SignalR&& out, const SignalT& in) -> decltype(impl::Fft(AsView(out), AsView(in))) {
	return impl::Fft(AsView(out), AsView(in));
}

template <class SignalR, class SignalT>
auto Ifft(SignalR&& out, const SignalT& in) -> decltype(impl::Ifft(AsView(out), AsView(in))) {
	return impl::Ifft(AsView(out), AsView(in));
}


template <class SignalT>
auto Fft(const SignalT& in, impl::FftFull) -> decltype(impl::Fft(AsView(in), FFT_FULL)) {
	return impl::Fft(AsView(in), FFT_FULL);
}

template <class SignalT>
auto Fft(const SignalT& in, impl::FftHalf) -> decltype(impl::Fft(AsView(in), FFT_HALF)) {
	return impl::Fft(AsView(in), FFT_HALF);
}

template <class SignalT>
auto Fft(const SignalT& in) -> decltype(impl::Fft(AsView(in))) {
	return impl::Fft(AsView(in));
}

template <class SignalT>
auto Ifft(const SignalT& in, impl::FftFull) -> decltype(impl::Ifft(AsView(in), FFT_FULL)) {
	return impl::Ifft(AsView(in), FFT_FULL);
}

template <class SignalT>
auto Ifft(const SignalT& in, impl::FftHalf, bool even) -> decltype(impl::Ifft(AsView(in), FFT_HALF, even)) {
	return impl::Ifft(AsView(in), FFT_HALF, even);
}

template <class SignalT>
auto Ifft(const SignalT& in) -> decltype(impl::Ifft(AsView(in))) {
	return impl::Ifft(AsView(in));
}


//------------------------------------------------------------------------------
// Utilities
//------------------------------------------------------------------------------

inline double FourierBin2Frequency(size_t binIdx, size_t numBins, uint64_t sampleRate) {
	return double(binIdx) / double(numBins) * double(sampleRate);
}

inline size_t FourierFrequency2Bin(double frequency, size_t numBins, uint64_t sampleRate) {
	return size_t(std::round(frequency / double(sampleRate) * double(numBins)));
}

namespace impl {
	template <class SignalR, class SignalT, std::enable_if_t<is_same_domain_v<SignalR, SignalT>, int> = 0>
	void BasicShift(SignalR&& out, const SignalT& in, size_t shift) {
		assert(out.Size() == in.Size());
		if (static_cast<const void*>(std::addressof(out)) == static_cast<const void*>(std::addressof(in))) {
			const auto first = out.begin();
			const auto mid = out.begin() + shift;
			const auto last = out.end();
			std::rotate(first, mid, last);
		}
		else {
			const auto first = in.begin();
			const auto mid = in.begin() + shift;
			const auto last = in.end();
			const auto write = out.begin();
			std::rotate_copy(first, mid, last, write);
		}
	}
} // namespace impl

template <class SignalR, class SignalT, std::enable_if_t<is_same_domain_v<SignalR, SignalT>, int> = 0>
void FftShift(SignalR&& out, const SignalT& in) {
	const size_t shift = (1 + out.Size()) / 2;
	impl::BasicShift(out, in, shift);
}

template <class SignalT, std::enable_if_t<is_signal_like_v<SignalT>, int> = 0>
SignalT FftShift(const SignalT& in) {
	SignalT out(in.Size());
	FftShift(out, in);
	return out;
}

template <class SignalR, class SignalT, std::enable_if_t<is_same_domain_v<SignalR, SignalT>, int> = 0>
void IfftShift(SignalR&& out, const SignalT& in) {
	const size_t shift = out.Size() / 2;
	impl::BasicShift(out, in, shift);
}

template <class SignalT, std::enable_if_t<is_signal_like_v<SignalT>, int> = 0>
SignalT IfftShift(const SignalT& in) {
	SignalT out(in.Size());
	IfftShift(out, in);
	return out;
}


} // namespace dspbb
