#pragma once

#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Primitives/SignalTraits.hpp>
#include <dspbb/Primitives/SignalView.hpp>
#include <dspbb/Filtering/FIR.hpp>


namespace dspbb {


template <class SignalR, class SignalT, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
void HilbertFirFromHalfbandIII(SignalR& out, const SignalT& halfband) {
	assert(halfband.Size() % 2 == 1);

	using R = typename std::decay_t<SignalR>::value_type;
	using T = typename std::decay_t<SignalT>::value_type;
	constexpr auto Domain = signal_traits<std::decay_t<SignalR>>::domain;

	const size_t filterSize = halfband.Size();
	constexpr size_t kernelSize = 32;
	static constexpr std::array<T, kernelSize> kernel = { 2, 0, -2, 0, 2, 0, -2, 0,
														  2, 0, -2, 0, 2, 0, -2, 0,
														  2, 0, -2, 0, 2, 0, -2, 0,
														  2, 0, -2, 0, 2, 0, -2, 0 };
	constexpr size_t kernelCenter = kernelSize / 2 - 1;
	constexpr size_t maxSizeSingleStep = kernelSize - 1;

	if (halfband.Size() <= maxSizeSingleStep) {
		const size_t offset = kernelCenter - filterSize / 2;
		SignalView<const T, Domain> kernelRegion{ kernel.data() + offset, kernel.data() + filterSize };
		Multiply(out, halfband, kernelRegion);
	}
	else {
		size_t tap = (filterSize / 2 - kernelCenter) % kernelSize;

		Multiply(SignalView<R, Domain>{ out.Data(), tap },
				 SignalView<const T, Domain>{ halfband.Data(), tap },
				 SignalView<const T, Domain>{ kernel.begin() + kernelSize - tap, kernel.end() });
		for (; tap < filterSize; tap += kernelSize) {
			Multiply(SignalView<R, Domain>{ out.Data() + tap, kernelSize },
					 SignalView<const T, Domain>{ halfband.Data() + tap, kernelSize },
					 SignalView<const T, Domain>{ kernel.data(), kernelSize });
		}
		const size_t lastChunkSize = filterSize - tap;
		Multiply(SignalView<R, Domain>{ out.Data() + tap, lastChunkSize },
				 SignalView<const T, Domain>{ halfband.Data() + tap, lastChunkSize },
				 SignalView<const T, Domain>{ kernel.data(), lastChunkSize });
	}
}

template <class SignalR, class SignalT, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void HilbertFirFromHalfbandIV(SignalR& out, const SignalT& halfband) {
	assert(halfband.Size() % 2 == 1);
}


template <class SignalR, class WindowFunc, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void HilbertFirWinIII(SignalR& out, WindowFunc windowFunc = windows::hamming) {
	FirLowpassWin(out, 0.5, windowFunc);
	HilbertFirFromHalfbandIII(out, out);
}

template <class SignalR, class SignalT, class WindowFunc, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void HilbertFirWinIII(SignalR& out, const SignalT& window) {
	FirLowpassWin(out, 0.5, window);
	HilbertFirFromHalfbandIII(out, out);
}

template <class T, eSignalDomain Domain, class WindowFunc>
auto HilbertFirWinIII(size_t taps, WindowFunc windowFunc = windows::hamming) {
	auto out = FirLowpassWin<T, Domain>(0.5, taps, windowFunc);
	HilbertFirFromHalfbandIII(out, out);
	return out;
}

template <class T, eSignalDomain Domain, class SignalT, class WindowFunc>
auto HilbertFirWinIII(size_t taps, const SignalT& window) {
	auto out = FirLowpassWin<T, Domain>(0.5, taps, window);
	HilbertFirFromHalfbandIII(out, out);
	return out;
}




} // namespace dspbb
