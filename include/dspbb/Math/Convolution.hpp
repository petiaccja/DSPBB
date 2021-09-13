#pragma once

#include "../ComputeKernels/Convolution.hpp"
#include "../Math/Convolution.hpp"
#include "../Math/DotProduct.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"

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
	const auto& shorter = mm.first;
	const auto& longer = mm.second;
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
	const auto& shorter = mm.first;
	const auto& longer = mm.second;
	return longer + shorter - 1;
}

namespace impl {
	template <class R, class T, class U, eSignalDomain Domain>
	auto ConvolutionOrdered(SignalView<R, Domain> r, SignalView<const T, Domain> u, SignalView<const U, Domain> v, size_t offset) {
		const size_t len = r.Length();
		Convolution(r.Data(), u.Data(), v.Data(), u.Size(), v.Size(), offset, len);
	}
} // namespace impl


// WARNING: offset is not tested for this, apart from full and central.
template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto Convolution(const SignalT& u, const SignalU& v, size_t offset, size_t length) {
	using T = typename signal_traits<std::decay_t<SignalT>>::type;
	using U = typename signal_traits<std::decay_t<SignalU>>::type;
	using R = decltype(std::declval<T>() * std::declval<U>());
	constexpr eSignalDomain Domain = signal_traits<std::decay_t<SignalT>>::domain;

	Signal<R, Domain> r(size_t(length), R(0));
	u.Size() >= v.Size() ? impl::ConvolutionOrdered(AsView(r), AsView(u), AsView(v), offset) :
							 impl::ConvolutionOrdered(AsView(r), AsView(v), AsView(u), offset);
	return r;
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto Convolution(const SignalT& u, const SignalU& v, convolution::impl::Full) {
	const size_t length = ConvolutionLength(u.Length(), v.Length(), convolution::full);
	size_t offset = 0;
	return Convolution(u, v, offset, length);
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto Convolution(const SignalT& u, const SignalU& v, convolution::impl::Central) {
	const size_t length = ConvolutionLength(u.Length(), v.Length(), convolution::central);
	size_t offset = std::min(u.Size() - 1, v.Size() - 1);
	return Convolution(u, v, offset, length);
}

} // namespace dspbb
