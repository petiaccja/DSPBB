#pragma once

#include "../../Primitives/SignalTraits.hpp"
#include "../../Utility/Numbers.hpp"

#include <array>
#include <cassert>
#include <cmath>


namespace dspbb::fir {


template <class SignalR, class SignalT, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
void MirrorResponse(SignalR&& mirrored, const SignalT& filter) {
	assert(mirrored.size() == filter.size());
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	T sign = T(1);
	for (size_t i = 0; i < filter.size(); ++i, sign *= T(-1)) {
		mirrored[i] = R(sign * filter[i]);
	}
}

template <class SignalR, class SignalT, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
void ComplementaryResponse(SignalR&& complementary, const SignalT& filter) {
	assert(filter.size() % 2 == 1);
	using R = typename signal_traits<std::decay_t<SignalR>>::type;
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	Multiply(complementary, filter, T(-1));
	complementary[complementary.size() / 2] += R(1);
}

template <class SignalR, class SignalT, class U, std::enable_if_t<is_mutable_signal_v<SignalR> && is_same_domain_v<SignalR, SignalT>, int> = 0>
void ShiftResponse(SignalR&& moved, const SignalT& filter, U normalizedFrequency) {
	assert(moved.size() == filter.size());
	const auto offset = static_cast<U>(filter.size() / 2);
	const U scale = pi_v<U> * normalizedFrequency;
	const size_t size = filter.size();
	for (size_t i = 0; i < size / 2; ++i) {
		const U x = (U(i) - offset) * scale;
		const U c = std::cos(x);
		moved[i] = c * filter[i];
		moved[size - i - 1] = c * filter[size - i - 1];
	}
	moved *= typename signal_traits<SignalT>::type(2);
}


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
void HalfbandToHilbertOdd(SignalR& out, const SignalT& halfband) {
	assert(halfband.size() % 2 == 1);
	assert(out.size() == halfband.size());

	using impl::kernelSize;
	using R = typename std::decay_t<SignalR>::value_type;
	using T = typename std::decay_t<SignalT>::value_type;
	constexpr auto Domain = signal_traits<std::decay_t<SignalR>>::domain;
	constexpr size_t kernelCenter = kernelSize / 2 - 1;
	constexpr size_t maxSizeSingleStep = kernelSize - 1;
	const BasicSignalView<const T, Domain> kernel(impl::kernel<T>.begin(), impl::kernel<T>.end());

	const size_t filterSize = halfband.size();

	if (halfband.size() <= maxSizeSingleStep) {
		const size_t offset = kernelCenter - filterSize / 2;
		const auto kernelRegion = kernel.subsignal(offset, filterSize);
		Multiply(out, halfband, kernelRegion);
	}
	else {
		size_t tap = (filterSize / 2 - kernelCenter) % kernelSize;

		Multiply(BasicSignalView<R, Domain>{ out.data(), tap },
				 BasicSignalView<const T, Domain>{ halfband.data(), tap },
				 kernel.subsignal(kernelSize - tap));
		for (; tap + kernelSize < filterSize; tap += kernelSize) {
			Multiply(BasicSignalView<R, Domain>{ out.data() + tap, kernelSize },
					 BasicSignalView<const T, Domain>{ halfband.data() + tap, kernelSize },
					 kernel);
		}
		const size_t lastChunkSize = filterSize - tap;
		Multiply(BasicSignalView<R, Domain>{ out.data() + tap, lastChunkSize },
				 BasicSignalView<const T, Domain>{ halfband.data() + tap, lastChunkSize },
				 kernel.subsignal(0, lastChunkSize));
	}
}

template <class SignalR, class SignalT, std::enable_if_t<is_mutable_signal_v<SignalR>, int> = 0>
void HalfbandToHilbertEven(SignalR& out, const SignalT& halfband) {
	assert(out.size() % 2 == 0);
	assert(out.size() * 2 - 1 == halfband.size());

	using impl::kernelSize;
	using R = typename std::decay_t<SignalR>::value_type;
	using T = typename std::decay_t<SignalT>::value_type;
	constexpr auto Domain = signal_traits<std::decay_t<SignalR>>::domain;
	constexpr size_t kernelCenter = kernelSize / 2 - 1;
	constexpr size_t maxSizeSingleStep = kernelSize - 1;

	std::array<T, kernelSize> scratchStorage;
	const BasicSignalView<T, Domain> scratch(scratchStorage.begin(), scratchStorage.end());
	const BasicSignalView<const T, Domain> kernel(impl::kernel<T>.begin(), impl::kernel<T>.end());

	const size_t filterSize = halfband.size();

	if (halfband.size() <= maxSizeSingleStep) {
		const size_t offset = kernelCenter - filterSize / 2;
		const auto kernelRegion = kernel.subsignal(offset, filterSize);
		const auto scratchRegion = scratch.subsignal(0, filterSize);
		Multiply(scratchRegion, halfband, kernelRegion);
		Decimate(out, scratchRegion, 2);
	}
	else {
		size_t tap = (filterSize / 2 - kernelCenter) % kernelSize;

		Multiply(scratch.subsignal(0, tap),
				 BasicSignalView<const T, Domain>{ halfband.data(), tap },
				 kernel.subsignal(kernelSize - tap));
		Decimate(BasicSignalView<T, Domain>(out.begin(), (tap + 1) / 2), scratch.subsignal(0, tap), 2);

		for (; tap + kernelSize < filterSize; tap += kernelSize) {
			Multiply(scratch,
					 BasicSignalView<const T, Domain>{ halfband.data() + tap, kernelSize },
					 kernel);
			Decimate(BasicSignalView<R, Domain>{ out.begin() + (tap + 1) / 2, (kernelSize + 1) / 2 }, scratch, 2);
		}

		const size_t lastChunkSize = filterSize - tap;
		Multiply(scratch.subsignal(0, lastChunkSize),
				 BasicSignalView<const T, Domain>{ halfband.data() + tap, lastChunkSize },
				 kernel.subsignal(0, lastChunkSize));
		Decimate(BasicSignalView<R, Domain>{ out.begin() + (tap + 1) / 2, (lastChunkSize + 1) / 2 }, scratch.subsignal(0, lastChunkSize), 2);
	}
}

} // namespace dspbb::fir