#pragma once

#include "../Math/DotProduct.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "../Utility/TypeTraits.hpp"

#include <complex>

namespace dspbb {

namespace convolution {
	namespace impl {
		class Central {};
		class Full {};
	} // namespace impl
	constexpr impl::Central central;
	constexpr impl::Full full;
} // namespace convolution


namespace impl {

	template <class T, class U>
	using ResultT = std::conditional_t<is_complex<T>::value || is_complex<U>::value,
									   std::complex<ProductT<remove_complex_t<T>, remove_complex_t<U>>>,
									   ProductT<remove_complex_t<T>, remove_complex_t<U>>>;

	template <class T, class U, eSignalDomain Domain>
	auto ConvolutionOrdered(SignalView<const T, Domain> base, SignalView<const U, Domain> running, convolution::impl::Full) {
		using R = ResultT<T, U>;

		assert(!base.Empty());
		assert(!running.Empty());
		assert(base.Size() > running.Size());

		const size_t paddingLength = running.Size() - 1;

		Signal<R, Domain> out;

		size_t outLength = base.Size() + paddingLength;
		out.Reserve(base.Size() + 2 * paddingLength);
		out.Append(Signal<R, Domain>(paddingLength, R(0)));
		out.Insert(out.end(), base.begin(), base.end());
		out.Append(Signal<R, Domain>(paddingLength, R(0)));
		// In-place convolution
		for (size_t offset = 0; offset < outLength; ++offset) {
			out[offset] = DotProduct(SignalView<const R, Domain>{ out.begin() + offset, out.end() }, running, running.Size());
		}
		out.Resize(outLength);

		return out;
	}

	template <class T, class U, eSignalDomain Domain>
	auto ConvolutionOrdered(SignalView<const T, Domain> base, SignalView<const U, Domain> running, convolution::impl::Central) {
		using R = ResultT<T, U>;

		assert(!base.Empty());
		assert(!running.Empty());
		assert(base.Size() > running.Size());

		const size_t paddingLength = running.Size() - 1;
		const size_t centralLength = base.Size() - paddingLength;

		Signal<R, Domain> out;

		out.Resize(centralLength);
		for (size_t offset = 0; offset < centralLength; ++offset) {
			out[offset] = DotProduct(base.SubSignal(offset), running, running.Size());
		}

		return out;
	}

} // namespace impl


/// <summary>
/// Ordinary convolution, EXCEPT it does not flip <paramref name="v"/> on the X axis.
/// </summary>
/// <typeparam name="T"> Any float or complex. </typeparam>
/// <typeparam name="U"> Any float or complex, operations should work with <typeparamref name="T"/>. </typeparam>
/// <typeparam name="PaddingMode"> One of <see cref="dspbb::convolution::full"/> or <see cref="dspbb::convolution::central"/>. </typeparam>
/// <param name="u"> The first argument of the convolution. </param>
/// <param name="v"> The second argument of the convolution. </param>
/// <returns> The result of the convolution. </returns>
template <class T, class U, eSignalDomain Domain, class PaddingMode>
auto ConvolutionFast(SignalView<const T, Domain> u, SignalView<const U, Domain> v, PaddingMode) {
	assert(!u.Empty());
	assert(!v.Empty());

	return u.Size() >= v.Size() ? impl::ConvolutionOrdered(u, v, PaddingMode{}) : impl::ConvolutionOrdered(v, u, PaddingMode{});
}

/// <summary>
/// Ordinary convolution, EXCEPT it does not flip <paramref name="v"/> on the X axis.
/// </summary>
/// <typeparam name="SignalT"> Either a Signal or SignalView. </typeparam>
/// <typeparam name="SignalU"> Either a Signal or SignalView, same domain as SignalT. </typeparam>
/// <typeparam name="PaddingMode"> One of <see cref="dspbb::convolution::full"/> or <see cref="dspbb::convolution::central"/>. </typeparam>
/// <param name="u"> The first argument of the convolution. </param>
/// <param name="v"> The second argument of the convolution. </param>
/// <returns> The result of the convolution. </returns>
template <class SignalT, class SignalU, class PaddingMode>
auto ConvolutionFast(const SignalT& u, const SignalU& v, PaddingMode) {
	return ConvolutionFast(AsConstView(u), AsConstView(v), PaddingMode{});
}


/// <summary> Calculates the length of the result of the convolution U*V. </summary>
/// <param name="lengthU"> Length of U. </param>
/// <param name="lengthV"> Length of V. </param>
inline size_t ConvolutionLength(size_t lengthU, size_t lengthV, convolution::impl::Central) {
	assert(lengthU > 0);
	assert(lengthV > 0);
	const auto& mm = std::minmax(lengthU, lengthV);
	const auto& shorter = mm.first, longer = mm.second;
	return longer - shorter + 1;
}

/// <summary> Calculates the length of the result of the convolution U*V. </summary>
/// <param name="lengthU"> Length of U. </param>
/// <param name="lengthV"> Length of V. </param>
inline size_t ConvolutionLength(size_t lengthU, size_t lengthV, convolution::impl::Full) {
	assert(lengthU > 0);
	assert(lengthV > 0);
	const auto& mm = std::minmax(lengthU, lengthV);
	const auto &shorter = mm.first, longer = mm.second;
	return longer + shorter - 1;
}


} // namespace dspbb
