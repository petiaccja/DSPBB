#pragma once

#include "../Math/Convolution.hpp"
#include "../Math/FFT.hpp"
#include "../Math/Solvers.hpp"
#include "../Utility/Interval.hpp"

#include <cmath>

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

		// Cost of doing OLA with fftSize=K, filterSize=F, and signal length N:
		// N/(K-F) * (2*k1*K log K + k2*K + k3*K)
		// Where k1, k2, and k3 are the constants for FFT, ADD, and MUL operations
		// from the big O notation. I.e. cost of FFT is O(K log K) = k1*K log K.
		// We need to equate the derivative of the cost function to zero to get the optimal FFT size.

		// Numerator of the derivative of the cost function. (d / d fftSize)
		inline double CostDX(double fftSize, double filterSize, double constFft, double constAdd, double constMul) {
			return filterSize * (2 * constFft + constAdd + constMul) + 2.0 * constFft * (filterSize * std::log(fftSize) - fftSize);
		}

		// Derivative of the above numerator.
		inline double CostD2X2(double fftSize, double filterSize, double constFft) {
			return 2.0 * constFft * (filterSize / fftSize - 1.0);
		}

		// No damn clue about these constants.
		// They depend on the CPU as well as the FFT algorithm and the vectorization of VMULPS and VADDPS.
		// Underestimating the constant for FFT and overestimating for MUL and ADD is less of an issue.
		// ^ That will suggest larger FFT than optimal, with a small performance hit.
		// Besides, if the user needs performance, he can benchmark the whole OLA and convolution himself!
		constexpr double kFft = 6.0;
		constexpr double kAdd = 1.0;
		constexpr double kMul = 3.0;

		// We can solve OlaCostDX = 0 with Newton's Method.
		inline double OptimalTheoreticalSize(double filterSize, double constFft = kFft, double constAdd = kAdd, double constMul = kMul) {
			auto myCostDX = [=](double fftSize) {
				return CostDX(fftSize, filterSize, constFft, constAdd, constMul);
			};
			auto myCostD2X2 = [=](double fftSize) {
				return CostD2X2(fftSize, filterSize, constFft);
			};

			const double x0 = 3.0 * filterSize; // d^2 / dx^2 is 0 at fftSize=filterSize, anything higher than filterSize is theoretically fine.
			return NewtonRaphson(myCostDX, myCostD2X2, x0);
		}

		inline size_t NextPowerOfTwo(size_t n) {
			assert(n < size_t(1) << (sizeof(size_t) * 8 - 1));
			if (n == 0) {
				return 0;
			}
			size_t p = 1;
			while (p < n) {
				p <<= 1;
			}
			return p;
		}

		inline size_t OptimalPracticalSize(size_t signalSize, size_t filterSize, double constFft = kFft, double constAdd = kAdd, double constMul = kMul) {
			size_t maxUsefulSize = ConvolutionLength(signalSize, filterSize, CONV_FULL);
			size_t suggestedSize = NextPowerOfTwo(size_t(OptimalTheoreticalSize(double(filterSize), constFft, constAdd, constMul)));
			if (suggestedSize * 3 / 4 < maxUsefulSize) {
				return suggestedSize;
			}
			return maxUsefulSize;
		}
	} // namespace ola

} // namespace impl


/// <summary> Convolve the two signals using FFTs. </summary>
/// <param name="out"> The convolved signal is written here. </param>
/// <param name="u"> The first argument to the convolution. </param>
/// <param name="v"> The second argument to the convolution. </param>
/// <param name="offset"> Controls the starting point of the output. </param>
/// <param name="chunkSize"> The size of the FFT to use. </param>
/// <param name="clearOut"> Set to false if the output buffer is already zeroed. </param>
/// <remarks> The subset of the full convolution given by [offset, offset + out.size()).
///		is written into out. </remarks>
template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT, same_domain_as_r<SignalR> SignalU>
void OverlapAdd(SignalR&& out, const SignalT& u, const SignalU& v, size_t offset, size_t chunkSize = 0, bool clearOut = true) {
	if (u.size() < v.size()) {
		return OverlapAdd(out, v, u, offset, chunkSize, clearOut);
	}
	if (chunkSize == 0) {
		chunkSize = impl::ola::OptimalPracticalSize(u.size(), v.size());
	}
	assert(chunkSize >= 2 * v.size() - 1);
	const size_t fullLength = ConvolutionLength(u.size(), v.size(), CONV_FULL);
	assert(offset + out.size() <= fullLength && "Result is outside of full convolution, thus contains some true zeros. I mean, it's ok, but you are probably doing it wrong.");
	if (clearOut) {
		using R = scalar_type_t<std::decay_t<SignalR>>;
		std::fill(out.begin(), out.end(), R(remove_complex_t<R>(0)));
	}

	using T = std::remove_cv_t<scalar_type_t<std::decay_t<SignalT>>>;
	using U = std::remove_cv_t<scalar_type_t<std::decay_t<SignalU>>>;
	constexpr eSignalDomain Domain = domain_v<std::decay_t<SignalT>>;
	constexpr auto is_complex_t = std::integral_constant<bool, is_complex_v<T>>{};
	constexpr auto is_complex_u = std::integral_constant<bool, is_complex_v<U>>{};

	BasicSignal<U, Domain> filter(chunkSize, U(0));
	std::copy(v.begin(), v.end(), filter.begin());
	const auto filterFd = impl::ola::FftFilter(filter, is_complex_t, is_complex_u);

	const Interval outExtent{ intptr_t(offset), intptr_t(offset + out.size()) };
	const Interval uExtent{ intptr_t(0), intptr_t(u.size()) };
	const Interval loopInterval = Intersection(uExtent, EncompassingUnion(outExtent, outExtent + intptr_t(1) - intptr_t(v.size())));

	BasicSignal<T, Domain> workingChunk(chunkSize, T(0));
	Interval uInterval = { loopInterval.first, loopInterval.first + intptr_t(v.size()) };
	Interval outInterval = { loopInterval.first, loopInterval.first + intptr_t(chunkSize) };
	for (; !IsDisjoint(outInterval, outExtent); uInterval += intptr_t(v.size()), outInterval += intptr_t(v.size())) {
		Interval uValidInterval = Intersection(uInterval, uExtent);
		const auto fillFirst = std::copy(u.begin() + uValidInterval.first, u.begin() + uValidInterval.last, workingChunk.begin());
		std::fill(fillFirst, workingChunk.end(), T(0));

		const auto workingChunkFd = impl::ola::FftChunk(workingChunk, is_complex_t, is_complex_u);
		const auto filteredChunkFd = workingChunkFd * filterFd;
		const auto filteredChunk = impl::ola::IfftChunk(filteredChunkFd, is_complex_t, is_complex_u, chunkSize);

		Interval outValidInterval = Intersection(outInterval, outExtent) - intptr_t(offset);
		Interval chunkValidInterval = Intersection(outInterval, outExtent) - uInterval.first;

		AsView(out).subsignal(outValidInterval.first, outValidInterval.size()) += AsView(filteredChunk).subsignal(chunkValidInterval.first, chunkValidInterval.size());
	}
}


/// <summary> Convolve the two signals using FFTs. </summary>
/// <param name="out"> The convolved signal is written here. </param>
/// <param name="u"> The first argument to the convolution. </param>
/// <param name="v"> The second argument to the convolution. </param>
/// <param name="chunkSize"> The size of the FFT to use. </param>
/// <param name="clearOut"> Set to false if the output buffer is already zeroed. </param>
template <mutable_signal_or_view_r SignalR, same_domain_as_r<SignalR> SignalT, same_domain_as_r<SignalR> SignalU, eConvolutionMethod Method>
auto OverlapAdd(SignalR&& out, const SignalT& u, const SignalU& v, std::integral_constant<eConvolutionMethod, Method> method, size_t chunkSize = 0, bool clearOut = true) {
	const auto length = ConvolutionLength(u.size(), v.size(), method);
	const size_t offset = ConvolutionOffset(u.size(), v.size(), method);

	assert(out.size() == length && "Use ConvolutionLength to calculate output size properly.");

	OverlapAdd(out, u, v, offset, chunkSize, clearOut);
}


/// <summary> Convolve the two signals using FFTs. </summary>
/// <param name="u"> The first argument to the convolution. </param>
/// <param name="v"> The second argument to the convolution. </param>
/// <param name="offset"> Controls the starting point of the output. </param>
/// <param name="length"> Controls the length of the output. </param>
/// <param name="chunkSize"> The size of the FFT to use. </param>
/// <returns> The subset of the full convolution given by [offset, offset + length). </returns>
template <signal_or_view SignalT, same_domain_as<SignalT> SignalU>
auto OverlapAdd(const SignalT& u, const SignalU& v, size_t offset, size_t length, size_t chunkSize = 0) {
	constexpr eSignalDomain Domain = domain_v<std::decay_t<SignalT>>;
	using T = scalar_type_t<std::decay_t<SignalT>>;
	using U = scalar_type_t<std::decay_t<SignalU>>;
	using R = multiplies_result_t<T, U>;

	BasicSignal<R, Domain> out(length, R(remove_complex_t<R>(0)));
	OverlapAdd(out, u, v, offset, chunkSize, false);
	return out;
}


/// <summary> Convolve the two signals using FFTs. </summary>
/// <param name="u"> The first argument to the convolution. </param>
/// <param name="v"> The second argument to the convolution. </param>
/// <param name="chunkSize"> The size of the FFT to use. </param>
/// <returns> The full or central convolution. </returns>
template <signal_or_view SignalT, same_domain_as<SignalT> SignalU, eConvolutionMethod Method>
auto OverlapAdd(const SignalT& u, const SignalU& v, std::integral_constant<eConvolutionMethod, Method> method, size_t chunkSize = 0) {
	const size_t length = ConvolutionLength(u.size(), v.size(), method);
	const size_t offset = ConvolutionOffset(u.size(), v.size(), method);
	return OverlapAdd(u, v, offset, length, chunkSize);
}

} // namespace dspbb