#pragma once

#include "../Math/Convolution.hpp"
#include "../Math/FFT.hpp"


namespace dspbb {

namespace impl {

	namespace ola {

		template <class SignalU>
		auto FftFilter(const SignalU& filter, std::false_type, std::false_type) {
			return Fft(filter, FFT_HALF);
		}
		template <class SignalU>
		auto FftFilter(const SignalU& filter, std::false_type, std::true_type) {
			return Fft(filter);
		}
		template <class SignalU>
		auto FftFilter(const SignalU& filter, std::true_type, std::false_type) {
			return Fft(filter, FFT_FULL);
		}
		template <class SignalU>
		auto FftFilter(const SignalU& filter, std::true_type, std::true_type) {
			return Fft(filter);
		}


		template <class SignalT, bool S, bool F>
		auto FftChunk(const SignalT& chunk, std::integral_constant<bool, S>, std::integral_constant<bool, F>) {
			// Same thing as filter, just reverse the order of the complexness parameters.
			return FftFilter(chunk, std::integral_constant<bool, F>{}, std::integral_constant<bool, S>{});
		}

		template <class SpectrumT>
		auto IfftChunk(const SpectrumT& fft, std::false_type, std::false_type, size_t fftSize) {
			return Ifft(fft, FFT_HALF, fftSize % 2 == 0);
		}
		template <class SpectrumT>
		auto IfftChunk(const SpectrumT& fft, std::false_type, std::true_type, size_t) {
			return Ifft(fft);
		}
		template <class SpectrumT>
		auto IfftChunk(const SpectrumT& fft, std::true_type, std::false_type, size_t) {
			return Ifft(fft);
		}
		template <class SpectrumT>
		auto IfftChunk(const SpectrumT& fft, std::true_type, std::true_type, size_t) {
			return Ifft(fft);
		}

	} // namespace ola

} // namespace impl


template <class SignalR, class SignalT, class SignalU, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT, SignalU>, int> = 0>
void OverlapAdd(SignalR&& out, const SignalT& u, const SignalU& v, size_t offset, size_t chunkSize) {
	if (u.Size() < v.Size()) {
		return OverlapAdd(out, v, u, offset, chunkSize);
	}

	if (chunkSize < v.Size()) {
		throw std::invalid_argument("Chunk size must be at least the size of the filter.");
	}
	const size_t fullLength = ConvolutionLength(u.Length(), v.Length(), CONV_FULL);
	if (offset + out.Size() > fullLength) {
		throw std::out_of_range("Result is outside of full convolution, thus contains some true zeros. I mean, it's ok, but you are probably doing it wrong.");
	}

	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	using U = typename signal_traits<std::decay_t<SignalU>>::type;
	constexpr eSignalDomain Domain = signal_traits<std::decay_t<SignalT>>::domain;
	constexpr auto is_complex_t = std::integral_constant<bool, is_complex_v<T>>{};
	constexpr auto is_complex_u = std::integral_constant<bool, is_complex_v<U>>{};

	const size_t stepSize = chunkSize - v.Size();
	Signal<U, Domain> filterChunk(chunkSize, U(0));
	std::copy(v.begin(), v.end(), filterChunk.begin());
	const auto filterChunkFd = impl::ola::FftFilter(filterChunk, is_complex_t, is_complex_u);

	Signal<T, Domain> workingChunk(chunkSize, T(0));

	for (size_t inIdx = 0; inIdx < u.Size(); inIdx += stepSize) {
		const auto inChunk = AsView(u).SubSignal(inIdx, std::min(stepSize, u.Size() - inIdx));
		std::copy(inChunk.begin(), inChunk.end(), workingChunk.begin());
		std::fill(workingChunk.begin() + inChunk.Size(), workingChunk.end(), T(0));
		auto workingChunkFd = impl::ola::FftChunk(workingChunk, is_complex_t, is_complex_u);
		workingChunkFd *= filterChunkFd;
		const auto filteredChunk = impl::ola::IfftChunk(workingChunkFd, is_complex_t, is_complex_u, chunkSize);

		const intptr_t outFirstRaw = intptr_t(inIdx) - intptr_t(offset);
		const intptr_t outLast = std::min(outFirstRaw + intptr_t(chunkSize), intptr_t(out.Size()));
		const intptr_t outFirst = std::max(intptr_t(0), outFirstRaw);
		auto outChunk = AsView(out).SubSignal(outFirst, outLast - outFirst);
		outChunk += AsConstView(filteredChunk).SubSignal(outFirst - outFirstRaw, outLast - outFirst);
	}
}

template <class SignalR, class SignalT, class SignalU, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT, SignalU>, int> = 0>
void OverlapAdd(SignalR&& out, const SignalT& u, const SignalU& v, impl::ConvFull, size_t chunkSize) {
	const size_t fullLength = ConvolutionLength(u.Length(), v.Length(), CONV_FULL);
	if (out.Size() != fullLength) {
		throw std::invalid_argument("Use ConvolutionLength to calculate output size properly.");
	}
	size_t offset = 0;
	OverlapAdd(out, u, v, offset, chunkSize);
}

template <class SignalR, class SignalT, class SignalU, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT, SignalU>, int> = 0>
void OverlapAdd(SignalR&& out, const SignalT& u, const SignalU& v, impl::ConvCentral, size_t chunkSize) {
	const size_t centralLength = ConvolutionLength(u.Length(), v.Length(), CONV_CENTRAL);
	if (out.Size() != centralLength) {
		throw std::invalid_argument("Use ConvolutionLength to calculate output size properly.");
	}
	size_t offset = std::min(u.Size() - 1, v.Size() - 1);
	OverlapAdd(out, u, v, offset, chunkSize);
}


template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto OverlapAdd(const SignalT& u, const SignalU& v, size_t offset, size_t length, size_t chunkSize) {
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	using U = typename signal_traits<std::decay_t<SignalU>>::type;
	using R = product_type_t<T, U>;
	constexpr eSignalDomain Domain = signal_traits<std::decay_t<SignalT>>::domain;

	Signal<R, Domain> out(length, R(0));
	OverlapAdd(out, u, v, offset, chunkSize);
	return out;
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto OverlapAdd(const SignalT& u, const SignalU& v, impl::ConvFull, size_t chunkSize) {
	const size_t length = ConvolutionLength(u.Length(), v.Length(), CONV_FULL);
	size_t offset = 0;
	return OverlapAdd(u, v, offset, length, chunkSize);
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto OverlapAdd(const SignalT& u, const SignalU& v, impl::ConvCentral, size_t chunkSize) {
	const size_t length = ConvolutionLength(u.Length(), v.Length(), CONV_CENTRAL);
	size_t offset = std::min(u.Size() - 1, v.Size() - 1);
	return OverlapAdd(u, v, offset, length, chunkSize);
}

} // namespace dspbb