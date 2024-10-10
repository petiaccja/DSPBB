#pragma once

#include "../../Signal/Signal.hpp"
#include "../../Signal/SignalView.hpp"
#include "../../Signal/Traits.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <numbers>


namespace dspbb::fir {


/// <summary> Mirror the frequency response of a FIR filter. </summary>
/// <param name="mirrored"> The mirrored filter is written here. </param>
/// <param name="filter"> The filter who's response to mirror. </param>
/// <remarks> For example, a low-pass filter at cutoff of 0.6 is turned into a
///		high-pass filter as cutoff 0.4. </remarks>
template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT>
void MirrorResponse(SignalR&& mirrored, const SignalT& filter) {
	assert(!IsAliasing(mirrored, filter) || IsFullyAliasing(mirrored, filter));
	assert(mirrored.size() == filter.size());

	using R = scalar_type_t<std::decay_t<SignalR>>;
	using T = scalar_type_t<std::decay_t<SignalT>>;

	T sign = T(1);
	for (size_t i = 0; i < filter.size(); ++i, sign *= T(-1)) {
		mirrored[i] = R(sign * filter[i]);
	}
}


/// <summary> Mirror the frequency response of a FIR filter in-place. </summary>
template <mutable_signal_or_view_r SignalR>
void MirrorResponse(SignalR&& filter) {
	MirrorResponse(filter, filter);
}


/// <summary> Find the complement of a FIR filter. </summary>
/// <param name="complementary"> The complementary filter is written here. </param>
/// <param name="filter"> The filter who's complement to find. </param>
/// <remarks> The complementary, when added to the original filter, will
///		form an all-pass filter. </remarks>
template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT>
void ComplementaryResponse(SignalR&& complementary, const SignalT& filter) {
	assert(!IsAliasing(complementary, filter) || IsFullyAliasing(complementary, filter));
	assert(filter.size() % 2 == 1);

	using R = scalar_type_t<std::decay_t<SignalR>>;
	using T = scalar_type_t<std::decay_t<SignalT>>;

	Multiply(complementary, filter, T(-1));
	complementary[complementary.size() / 2] += R(1);
}


/// <summary> Find the complement of a FIR filter in-place. </summary>
template <mutable_signal_or_view_r SignalR>
void ComplementaryResponse(SignalR&& filter) {
	ComplementaryResponse(filter, filter);
}


/// <summary> Shift the response of a filter by a given frequency. </summary>
/// <param name="shifted"> The shifted filter is written here. </param>
/// <param name="filter"> The filter who's response to shift. </param>
/// <param name="normalizedFrequency"> The amount of shift in normalized frequency [-1, 1]. </param>
template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT, class U>
void ShiftResponse(SignalR&& shifted, const SignalT& filter, U normalizedFrequency) {
	assert(!IsAliasing(shifted, filter) || IsFullyAliasing(shifted, filter));
	assert(shifted.size() == shifted.size());

	const auto offset = static_cast<U>(filter.size() / 2);
	const U scale = std::numbers::pi_v<U> * normalizedFrequency;
	const size_t size = filter.size();
	for (size_t i = 0; i < size / 2; ++i) {
		const U x = (U(i) - offset) * scale;
		const U c = std::cos(x);
		shifted[i] = c * filter[i];
		shifted[size - i - 1] = c * filter[size - i - 1];
	}
	shifted *= scalar_type_t<SignalT>(2);
}


/// <summary> Shift the response of a filter by a given frequency in-place. </summary>
template <mutable_signal_or_view_r SignalR, class U>
void ShiftResponse(SignalR&& filter, U normalizedFrequency) {
	ShiftResponse(filter, filter, normalizedFrequency);
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


/// <summary> Convert a halfband (low-pass with cutoff at 0.5) filter to a hilbert filter. </summary>
/// <param name="out"> The hilbert filter is written here. </param>
/// <param name="halfband"> The halfband filter, which must have an odd number of coefficients. </param>
template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT>
void HalfbandToHilbertOdd(SignalR&& out, const SignalT& halfband) {
	assert(!IsAliasing(out, halfband) || IsFullyAliasing(out, halfband));
	assert(halfband.size() % 2 == 1);
	assert(out.size() == halfband.size());

	using impl::kernelSize;
	using R = typename std::decay_t<SignalR>::value_type;
	using T = typename std::decay_t<SignalT>::value_type;
	constexpr auto Domain = domain_v<std::decay_t<SignalR>>;
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


/// <summary> Convert a halfband (low-pass with cutoff at 0.5) filter to a hilbert filter in-place. </summary>
template <mutable_signal_or_view_r SignalR>
void HalfbandToHilbertOdd(SignalR&& filter) {
	HalfbandToHilbertOdd(filter, filter);
}


/// <summary> Convert a halfband (low-pass with cutoff at 0.5) filter to a hilbert filter. </summary>
/// <param name="out"> The hilbert filter is written here. The size must be len(halfband + 1) / 2. </param>
/// <param name="halfband"> The halfband filter, which must have an odd number of coefficients. </param>
template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT>
void HalfbandToHilbertEven(SignalR&& out, const SignalT& halfband) {
	assert(!IsAliasing(out, halfband));
	assert(out.size() % 2 == 0);
	assert(out.size() * 2 - 1 == halfband.size());

	using impl::kernelSize;
	using R = typename std::decay_t<SignalR>::value_type;
	using T = typename std::decay_t<SignalT>::value_type;
	constexpr auto Domain = domain_v<std::decay_t<SignalR>>;
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