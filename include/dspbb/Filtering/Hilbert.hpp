#pragma once

#include <dspbb/Filtering/FIR.hpp>
#include <dspbb/Filtering/Interpolation.hpp>
#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Primitives/SignalTraits.hpp>
#include <dspbb/Primitives/SignalView.hpp>


namespace dspbb {


namespace impl {

	constexpr size_t kernelSize = 32;

	template <class T>
	constexpr std::array<T, kernelSize> kernel = {
		2, 0, -2, 0, 2, 0, -2, 0,
		2, 0, -2, 0, 2, 0, -2, 0,
		2, 0, -2, 0, 2, 0, -2, 0,
		2, 0, -2, 0, 2, 0, -2, 0
	};
} // namespace impl


template <class SignalR, class SignalT, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
void HilbertFirFromHalfbandIII(SignalR& out, const SignalT& halfband) {
	assert(halfband.Size() % 2 == 1);

	using impl::kernelSize;
	using R = typename std::decay_t<SignalR>::value_type;
	using T = typename std::decay_t<SignalT>::value_type;
	constexpr auto Domain = signal_traits<std::decay_t<SignalR>>::domain;
	constexpr size_t kernelCenter = kernelSize / 2 - 1;
	constexpr size_t maxSizeSingleStep = kernelSize - 1;
	const SignalView<const T, Domain> kernel(impl::kernel<T>.begin(), impl::kernel<T>.end());

	const size_t filterSize = halfband.Size();

	if (halfband.Size() <= maxSizeSingleStep) {
		const size_t offset = kernelCenter - filterSize / 2;
		const auto kernelRegion = kernel.SubSignal(offset, filterSize);
		Multiply(out, halfband, kernelRegion);
	}
	else {
		size_t tap = (filterSize / 2 - kernelCenter) % kernelSize;

		Multiply(SignalView<R, Domain>{ out.Data(), tap },
				 SignalView<const T, Domain>{ halfband.Data(), tap },
				 kernel.SubSignal(kernelSize - tap));
		for (; tap + kernelSize < filterSize; tap += kernelSize) {
			Multiply(SignalView<R, Domain>{ out.Data() + tap, kernelSize },
					 SignalView<const T, Domain>{ halfband.Data() + tap, kernelSize },
					 kernel);
		}
		const size_t lastChunkSize = filterSize - tap;
		Multiply(SignalView<R, Domain>{ out.Data() + tap, lastChunkSize },
				 SignalView<const T, Domain>{ halfband.Data() + tap, lastChunkSize },
				 kernel.SubSignal(0, lastChunkSize));
	}
}

template <class SignalR, class SignalT, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void HilbertFirFromHalfbandIV(SignalR& out, const SignalT& halfband) {
	assert(out.Size() * 2 - 1 == halfband.Size());
	assert(out.Size() % 2 == 0);

	using impl::kernelSize;
	using R = typename std::decay_t<SignalR>::value_type;
	using T = typename std::decay_t<SignalT>::value_type;
	constexpr auto Domain = signal_traits<std::decay_t<SignalR>>::domain;
	constexpr size_t kernelCenter = kernelSize / 2 - 1;
	constexpr size_t maxSizeSingleStep = kernelSize - 1;

	std::array<T, kernelSize> scratchStorage;
	const SignalView<T, Domain> scratch(scratchStorage.begin(), scratchStorage.end());
	const SignalView<const T, Domain> kernel(impl::kernel<T>.begin(), impl::kernel<T>.end());

	const size_t filterSize = halfband.Size();

	if (halfband.Size() <= maxSizeSingleStep) {
		const size_t offset = kernelCenter - filterSize / 2;
		const auto kernelRegion = kernel.SubSignal(offset, filterSize);
		const auto scratchRegion = scratch.SubSignal(0, filterSize);
		Multiply(scratchRegion, halfband, kernelRegion);
		Decimate(out, scratchRegion, 2);
	}
	else {
		size_t tap = (filterSize / 2 - kernelCenter) % kernelSize;

		Multiply(scratch.SubSignal(0, tap),
				 SignalView<const T, Domain>{ halfband.Data(), tap },
				 kernel.SubSignal(kernelSize - tap));
		Decimate(SignalView<T, Domain>(out.begin(), (tap + 1) / 2), scratch.SubSignal(0, tap), 2);

		for (; tap + kernelSize < filterSize; tap += kernelSize) {
			Multiply(scratch,
					 SignalView<const T, Domain>{ halfband.Data() + tap, kernelSize },
					 kernel);
			Decimate(SignalView<R, Domain>{ out.begin() + (tap + 1) / 2, (kernelSize + 1) / 2 }, scratch, 2);
		}

		const size_t lastChunkSize = filterSize - tap;
		Multiply(scratch.SubSignal(0, lastChunkSize),
				 SignalView<const T, Domain>{ halfband.Data() + tap, lastChunkSize },
				 kernel.SubSignal(0, lastChunkSize));
		Decimate(SignalView<R, Domain>{ out.begin() + (tap + 1) / 2, (lastChunkSize + 1) / 2 }, scratch.SubSignal(0, lastChunkSize), 2);
	}
}


template <class SignalR, class WindowFunc, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void HilbertFirWinIII(SignalR& out, WindowFunc windowFunc = windows::hamming) {
	FirFilter(out, Lowpass(0.5), Windowed(windowFunc));
	HilbertFirFromHalfbandIII(out, out);
}

template <class SignalR, class SignalT, class WindowFunc, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void HilbertFirWinIII(SignalR& out, const SignalT& window) {
	FirFilter(out, Lowpass(0.5), Windowed(window));
	HilbertFirFromHalfbandIII(out, out);
}

template <class T, eSignalDomain Domain, class WindowFunc>
auto HilbertFirWinIII(size_t taps, WindowFunc windowFunc = windows::hamming) {
	auto out = FirFilter<T, Domain>(taps, Lowpass(0.5), Windowed(windowFunc));
	HilbertFirFromHalfbandIII(out, out);
	return out;
}

template <class T, eSignalDomain Domain, class SignalT, class WindowFunc>
auto HilbertFirWinIII(size_t taps, const SignalT& window) {
	auto out = FirFilter<T, Domain>(taps, Lowpass(0.5), Windowed(window));
	HilbertFirFromHalfbandIII(out, out);
	return out;
}



template <class SignalR, class WindowFunc, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void HilbertFirWinIV(SignalR& out, WindowFunc windowFunc = windows::hamming) {
	FirFilter(out, Lowpass(0.5), Windowed(windowFunc));
	HilbertFirFromHalfbandIV(out, out);
}

template <class SignalR, class SignalT, class WindowFunc, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void HilbertFirWinIV(SignalR& out, const SignalT& window) {
	FirFilter(out, Lowpass(0.5), Windowed(window));
	HilbertFirFromHalfbandIV(out, out);
}

template <class T, eSignalDomain Domain, class WindowFunc>
auto HilbertFirWinIV(size_t taps, WindowFunc windowFunc = windows::hamming) {
	const auto halfband = FirFilter<T, Domain>(2 * taps - 1, Lowpass(0.5), Windowed(windowFunc));
	Signal<T, Domain> out(taps);
	HilbertFirFromHalfbandIV(out, halfband);
	return out;
}

template <class T, eSignalDomain Domain, class SignalT, class WindowFunc>
auto HilbertFirWinIV(size_t taps, const SignalT& window) {
	const auto halfband = FirFilter<T, Domain>(2 * taps - 1, Lowpass(0.5), Windowed(window));
	Signal<T, Domain> out(taps);
	HilbertFirFromHalfbandIV(out, halfband);
	return out;
}



} // namespace dspbb
