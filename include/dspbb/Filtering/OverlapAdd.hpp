#pragma once

#include "Convolution.hpp"
#include "FFT.hpp"


namespace dspbb {

namespace impl {

	namespace ola {

		template <class SignalU>
		auto MakeFrequencyFilter(const SignalU& filter, std::false_type, std::false_type) {
			return FourierTransform(filter, false);
		}
		template <class SignalU>
		auto MakeFrequencyFilter(const SignalU& filter, std::false_type, std::true_type) {
			return FourierTransform(filter);
		}
		template <class SignalU>
		auto MakeFrequencyFilter(const SignalU& filter, std::true_type, std::false_type) {
			return FourierTransform(filter, true);
		}
		template <class SignalU>
		auto MakeFrequencyFilter(const SignalU& filter, std::true_type, std::true_type) {
			return FourierTransform(filter);
		}


		template <class SignalT, bool S, bool F>
		auto MakeFrequencyChunk(const SignalT& chunk, std::integral_constant<bool, S>, std::integral_constant<bool, F>) {
			// Same thing as filter, just reverse the order of the complexness parameters.
			return MakeFrequencyFilter(chunk, std::integral_constant<bool, F>{}, std::integral_constant<bool, S>{});
		}

		template <class SpectrumT>
		auto InvertChunk(const SpectrumT& fft, std::false_type, std::false_type, size_t fftSize) {
			return InverseFourierTransformR(fft, fftSize);
		}
		template <class SpectrumT>
		auto InvertChunk(const SpectrumT& fft, std::false_type, std::true_type, size_t fftSize) {
			return InverseFourierTransformC(fft);
		}
		template <class SpectrumT>
		auto InvertChunk(const SpectrumT& fft, std::true_type, std::false_type, size_t fftSize) {
			return InverseFourierTransformC(fft);
		}
		template <class SpectrumT>
		auto InvertChunk(const SpectrumT& fft, std::true_type, std::true_type, size_t fftSize) {
			return InverseFourierTransformC(fft);
		}

	} // namespace ola

} // namespace impl



template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto OverlapAdd(const SignalT& signal, const SignalU& filter, size_t chunkSize, size_t overlapSize, size_t offset, size_t length) {
	assert(chunkSize + overlapSize >= filter.Size());
	const size_t fftSize = chunkSize + overlapSize;

	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	using U = typename signal_traits<std::decay_t<SignalU>>::type;
	using R = decltype(std::declval<T>() * std::declval<U>());
	constexpr eSignalDomain Domain = signal_traits<std::decay_t<SignalT>>::domain;
	constexpr auto is_complex_t = std::integral_constant<bool, is_complex_v<T>>{};
	constexpr auto is_complex_u = std::integral_constant<bool, is_complex_v<U>>{};

	Signal<U, Domain> paddedFilter(fftSize, U(0));
	std::copy(filter.begin(), filter.end(), paddedFilter.begin());
	const auto frequencyFilter = impl::ola::MakeFrequencyFilter(paddedFilter, is_complex_t, is_complex_u);

	Signal<R, Domain> r(length, R(0));
	Signal<T, Domain> paddedChunk(fftSize, T(0));

	const auto signalView = AsConstView(signal);
	const auto rView = AsView(r);
	for (size_t i = 0; i < signal.Size(); i += chunkSize) {
		const auto chunkView = signalView.SubSignal(i, std::min(chunkSize, signalView.Size() - i));
		std::copy(chunkView.begin(), chunkView.end(), paddedChunk.begin());
		std::fill(paddedChunk.begin() + chunkView.Size(), paddedChunk.end(), T(0));
		auto frequencyChunk = impl::ola::MakeFrequencyChunk(paddedChunk, is_complex_t, is_complex_u);
		frequencyChunk *= frequencyFilter;
		const auto filteredChunk = impl::ola::InvertChunk(frequencyChunk, is_complex_t, is_complex_u, fftSize);
		
		const intptr_t outFirstRaw = intptr_t(i) - intptr_t(offset);
		const intptr_t outLast = std::min(outFirstRaw + intptr_t(fftSize), intptr_t(rView.Size()));
		const intptr_t outFirst = std::max(intptr_t(0), outFirstRaw);
		auto rChunkView = rView.SubSignal(outFirst, outLast - outFirst);
		rChunkView += AsConstView(filteredChunk).SubSignal(outFirst - outFirstRaw, outLast - outFirst);
	}

	return r;
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto OverlapAdd(const SignalT& u, const SignalU& v, size_t chunkSize, size_t overlapSize, convolution::impl::Full) {
	const size_t length = ConvolutionLength(u.Length(), v.Length(), convolution::full);
	size_t offset = 0;
	return OverlapAdd(u, v, chunkSize, overlapSize, offset, length);
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto OverlapAdd(const SignalT& u, const SignalU& v, size_t chunkSize, size_t overlapSize, convolution::impl::Central) {
	const size_t length = ConvolutionLength(u.Length(), v.Length(), convolution::central);
	size_t offset = std::min(u.Size() - 1, v.Size() - 1);
	return OverlapAdd(u, v, chunkSize, overlapSize, offset, length);
}

} // namespace dspbb