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


/// <summary> Calculates the length of the result of the convolution U*V. </summary>
/// <param name="lengthU"> Length of U. </param>
/// <param name="lengthV"> Length of V. </param>
inline size_t ConvolutionLength(size_t lengthU, size_t lengthV, convolution::impl::Central) {
	if (lengthU == 0 || lengthV == 0) {
		return 0;
	}
	const auto& mm = std::minmax(lengthU, lengthV);
	const auto &shorter = mm.first, longer = mm.second;
	return longer - shorter + 1;
}

/// <summary> Calculates the length of the result of the convolution U*V. </summary>
/// <param name="lengthU"> Length of U. </param>
/// <param name="lengthV"> Length of V. </param>
inline size_t ConvolutionLength(size_t lengthU, size_t lengthV, convolution::impl::Full) {
	if (lengthU == 0 || lengthV == 0) {
		return 0;
	}
	const auto& mm = std::minmax(lengthU, lengthV);
	const auto &shorter = mm.first, longer = mm.second;
	return longer + shorter - 1;
}

namespace impl {

	template <class T, class U>
	using ResultT = std::conditional_t<is_complex<T>::value || is_complex<U>::value,
									   std::complex<ProductT<remove_complex_t<T>, remove_complex_t<U>>>,
									   ProductT<remove_complex_t<T>, remove_complex_t<U>>>;

	template <class T, class U, eSignalDomain Domain>
	auto ConvolutionOrdered(SignalView<const T, Domain> u, SignalView<const U, Domain> v, convolution::impl::Full) {

		const intptr_t lenfull = (intptr_t)ConvolutionLength(u.Length(), v.Length(), convolution::full);
		const intptr_t lencentral = (intptr_t)ConvolutionLength(u.Length(), v.Length(), convolution::central);
		
		const intptr_t lenu = u.Length();
		const intptr_t lenv = v.Length();
		const intptr_t padding = lenv - 1;

		assert(lenu > lenv);

		using R = ResultT<T, U>;
		Signal<R, Domain> out(size_t(lenfull), R(0));

		for (intptr_t i = 0; i < padding; ++i) {
			for (intptr_t k = 0; k <= i; ++k) {
				out[i] += R(u[k]) * R(v[i - k]);
			}
		}
		for (intptr_t i = 0; i < lencentral; ++i) {
			for (intptr_t k = 0; k < lenv; ++k) {
				out[i + padding] += R(u[i + k]) * R(v[lenv - k - 1]);
			}
		}
		for (intptr_t i = 1; i < lenv; ++i) {
			for (intptr_t k = i; k < lenv; ++k) {
				out[i + padding + lencentral - 1] += R(u[lencentral + k - 1]) * R(v[lenv - k + i - 1]);
			}
		}

		return out;
	}

	template <class T, class U, eSignalDomain Domain>
	auto ConvolutionOrdered(SignalView<const T, Domain> u, SignalView<const U, Domain> v, convolution::impl::Central) {
		const intptr_t lenout = (intptr_t)ConvolutionLength(u.Length(), v.Length(), convolution::central);
		const intptr_t lenu = (intptr_t)u.Length();
		const intptr_t lenv = (intptr_t)v.Length();

		assert(lenu > lenv);

		using R = ResultT<T, U>;
		Signal<R, Domain> out(size_t(lenout), R(0));

		for (intptr_t k = 0; k < lenv; ++k) {
			for (intptr_t i = 0; i < lenout; ++i) {
				out[i] += R(u[i + k]) * R(v[lenv - k - 1]);
			}
		}

		return out;
	}

} // namespace impl


/// <summary>
/// Ordinary convolution.
/// </summary>
/// <typeparam name="T"> Any float or complex. </typeparam>
/// <typeparam name="U"> Any float or complex, operations should work with <typeparamref name="T"/>. </typeparam>
/// <typeparam name="PaddingMode"> One of <see cref="dspbb::convolution::full"/> or <see cref="dspbb::convolution::central"/>. </typeparam>
/// <param name="u"> The first argument of the convolution. </param>
/// <param name="v"> The second argument of the convolution. </param>
/// <returns> The result of the convolution. </returns>
template <class T, class U, eSignalDomain Domain, class PaddingMode>
auto Convolution(SignalView<const T, Domain> u, SignalView<const U, Domain> v, PaddingMode) {
	return u.Size() >= v.Size() ? impl::ConvolutionOrdered(u, v, PaddingMode{}) : impl::ConvolutionOrdered(v, u, PaddingMode{});
}

/// <summary>
/// Ordinary convolution.
/// </summary>
/// <typeparam name="SignalT"> Either a Signal or SignalView. </typeparam>
/// <typeparam name="SignalU"> Either a Signal or SignalView, same domain as SignalT. </typeparam>
/// <typeparam name="PaddingMode"> One of <see cref="dspbb::convolution::full"/> or <see cref="dspbb::convolution::central"/>. </typeparam>
/// <param name="u"> The first argument of the convolution. </param>
/// <param name="v"> The second argument of the convolution. </param>
/// <returns> The result of the convolution. </returns>
template <class SignalT, class SignalU, class PaddingMode>
auto Convolution(const SignalT& u, const SignalU& v, PaddingMode) {
	return Convolution(AsConstView(u), AsConstView(v), PaddingMode{});
}

} // namespace dspbb
