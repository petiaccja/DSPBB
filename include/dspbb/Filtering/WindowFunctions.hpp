#pragma once

#include "../Math/Statistics.hpp"
#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "../Utility/Algorithm.hpp"
#include "../Utility/Numbers.hpp"
#include "../Utility/TypeTraits.hpp"

#include <cmath>
#include <functional>


namespace dspbb {

//------------------------------------------------------------------------------
// Assess properties of windows.
//------------------------------------------------------------------------------
template <class T, eSignalDomain Domain>
T CoherentGain(SignalView<T, Domain> window) {
	return Sum(window) / remove_complex_t<T>(window.Size());
}

template <class T, eSignalDomain Domain>
T CoherentGain(const Signal<T, Domain>& window) {
	return CoherentGain(AsConstView(window));
}

template <class T, eSignalDomain Domain>
T EnergyGain(SignalView<T, Domain> window) {
	return SumSquare(window) / T(window.Size());
}

template <class T, eSignalDomain Domain>
T EnergyGain(const Signal<T, Domain>& window) {
	return EnergyGain(AsConstView(window));
}


//------------------------------------------------------------------------------
// List of window functions.
//------------------------------------------------------------------------------
template <class T, eSignalDomain Domain = eSignalDomain::TIME, std::enable_if_t<!std::is_const<T>::value, int> = 0>
void HammingWindow(SignalView<T, Domain> out) {
	using U = remove_complex_t<T>;

	const U N = U(out.Length());
	const U preSize = U(2) * pi_v<U> / (N - U(1));

	std::iota(out.begin(), out.end(), U(0.0));
	out *= preSize;
	Cos(out, out);
	out *= U(-0.46);
	out += U(0.54);
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME, std::enable_if_t<!std::is_const<T>::value, int> = 0>
void FlatTopWindow(SignalView<T, Domain> out) {
	using U = remove_complex_t<T>;

	U c0 = U(0.21557895);
	U c1 = U(-0.41663158);
	U c2 = U(0.277263158);
	U c3 = U(-0.083578947);
	U c4 = U(0.006947368);

	const U N = U(out.Length());
	U preSize1 = U(2) * pi_v<U> / (N - U(1));
	U preSize2 = U(4) * pi_v<U> / (N - U(1));
	U preSize3 = U(6) * pi_v<U> / (N - U(1));
	U preSize4 = U(8) * pi_v<U> / (N - U(1));

	std::iota(out.begin(), out.end(), U(0.0));
	Apply(out, [&](T k) -> T {
		U kreal = *reinterpret_cast<const U*>(&k); // We can do that safely for std::complex, it's a POD type.
		return c0
			   + c1 * std::cos(preSize1 * kreal)
			   + c2 * std::cos(preSize2 * kreal)
			   + c3 * std::cos(preSize3 * kreal)
			   + c4 * std::cos(preSize4 * kreal);
	});
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME, std::enable_if_t<!std::is_const<T>::value, int> = 0>
void RectangularWindow(SignalView<T, Domain> out) {
	std::fill(out.begin(), out.end(), T(1.0));
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
Signal<T, Domain> HammingWindow(size_t length) {
	Signal<T, Domain> window(length);
	HammingWindow(AsView(window));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
Signal<T, Domain> FlatTopWindow(size_t length) {
	Signal<T, Domain> window(length);
	FlatTopWindow(AsView(window));
	return window;
}

template <class T, eSignalDomain Domain = eSignalDomain::TIME>
Signal<T, Domain> RectangularWindow(size_t length) {
	Signal<T, Domain> window(length, T(1.0));
	return window;
}

//------------------------------------------------------------------------------
// Helper for when you have to pass a window function as an argument.
//------------------------------------------------------------------------------
namespace windows {
	namespace impl {
		template <class T, eSignalDomain Domain>
		using FunctionPtr = Signal<T, Domain> (*)(size_t);

		template <class T, eSignalDomain Domain>
		using Functional = std::function<Signal<T, Domain>(size_t)>;

		template <class T, eSignalDomain Domain>
		using FunctionPtrInplace = void (*)(SignalView<T, Domain>);

		template <class T, eSignalDomain Domain>
		using FunctionalInplace = std::function<void(SignalView<T, Domain>)>;
	} // namespace impl

#define DSPBB_WINDOW_CONVERSION_HELPER(VAR_NAME, FUNC_NAME)                                                         \
	struct FUNC_NAME##Conv {                                                                                        \
		template <class T, eSignalDomain Domain>                                                                    \
		operator impl::FunctionPtr<T, Domain>() {                                                                   \
			return static_cast<impl::FunctionPtr<T, Domain>>(FUNC_NAME);                                            \
		}                                                                                                           \
		template <class T, eSignalDomain Domain>                                                                    \
		operator impl::Functional<T, Domain>() {                                                                    \
			return impl::Functional<T, Domain>(static_cast<impl::FunctionPtr<T, Domain>>(FUNC_NAME));               \
		}                                                                                                           \
		template <class T, eSignalDomain Domain>                                                                    \
		operator impl::FunctionPtrInplace<T, Domain>() {                                                            \
			return static_cast<impl::FunctionPtrInplace<T, Domain>>(FUNC_NAME);                                     \
		}                                                                                                           \
		template <class T, eSignalDomain Domain>                                                                    \
		operator impl::FunctionalInplace<T, Domain>() {                                                             \
			return impl::FunctionalInplace<T, Domain>(static_cast<impl::FunctionPtrInplace<T, Domain>>(FUNC_NAME)); \
		}                                                                                                           \
	} static VAR_NAME;

	DSPBB_WINDOW_CONVERSION_HELPER(hamming, HammingWindow)
	DSPBB_WINDOW_CONVERSION_HELPER(flattop, FlatTopWindow)
	DSPBB_WINDOW_CONVERSION_HELPER(rectangular, RectangularWindow)
} // namespace windows



// TODO: do something about this
#if (defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || __cplusplus >= 201703L
template <class T, eSignalDomain Domain = eSignalDomain::TIME>
Signal<T, Domain> KaiserWindow(size_t length, T alpha = T(0.5)) {
	using UnderlyingT = remove_complex_t<T>;

	Signal<T, Domain> window;
	window.Resize(length);
	for (size_t i = 0; i < length; ++i) {
		constexpr UnderlyingT pi = UnderlyingT(3.141592653589793238);
		const UnderlyingT x = 2 * UnderlyingT(i) / UnderlyingT(length) - 1;
		const UnderlyingT value = UnderlyingT(std::cyl_bessel_i(UnderlyingT(0), pi * alpha * std::sqrt(1 - x * x)) / std::cyl_bessel_i(UnderlyingT(0), pi * alpha));
		window[i] = value;
	}
	return window;
}
#endif

} // namespace dspbb