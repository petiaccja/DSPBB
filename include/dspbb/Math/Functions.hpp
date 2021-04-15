#pragma once

#include <complex>
#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Primitives/SignalTraits.hpp>
#include <dspbb/Primitives/SignalView.hpp>
#include <dspbb/Utility/Algorithm.hpp>
#include <dspbb/Vectorization/ComplexFunctions.hpp>
#include <dspbb/Vectorization/MathFunctions.hpp>
#include <type_traits>


namespace dspbb {


#define DSPBB_IMPL_FUNCTION_2_PARAM(NAME, FUNC)                                                                                                                        \
	template <class SignalT, class SignalU, std::enable_if_t<is_mutable_signal_v<SignalT> && is_same_domain_v<std::decay_t<SignalT>, std::decay_t<SignalU>>, int> = 0> \
	auto NAME(SignalT&& out, const SignalU& in) {                                                                                                                      \
		return UnaryOperationVectorized(out.Data(), in.Data(), out.Length(), [](auto v) { return math_functions::FUNC(v); });                                          \
	}

#define DSPBB_IMPL_FUNCTION_1_PARAM(NAME)                                                        \
	template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0> \
	auto NAME(const SignalT& signal) {                                                           \
		SignalT r(signal.Size());                                                                \
		NAME(r, signal);                                                                         \
		return r;                                                                                \
	}

#define DSPBB_IMPL_FUNCTION(NAME, FUNC)     \
	DSPBB_IMPL_FUNCTION_2_PARAM(NAME, FUNC) \
	DSPBB_IMPL_FUNCTION_1_PARAM(NAME)


//------------------------------------------------------------------------------
// Complex number functions
//------------------------------------------------------------------------------

#define DSPBB_IMPL_COMPLEX_FUNCTION_2_PARAM(NAME, VECOP, OP, FUNC)                                                                                                     \
	template <class SignalT,                                                                                                                                           \
			  class SignalU,                                                                                                                                           \
			  class T,                                                                                                                                                 \
			  std::enable_if_t<is_mutable_signal_v<SignalT> && is_same_domain_v<std::decay_t<SignalT>, std::decay_t<SignalU>>, int> = 0>                               \
	auto NAME(SignalT&& out, const SignalU& in, int, std::complex<T>) {                                                                                                \
                                                                                                                                                                       \
		return UnaryOperationVectorized(out.Data(),                                                                                                                    \
										in.Data(),                                                                                                                     \
										out.Length(),                                                                                                                  \
										complex_functions::VECOP<T>::stride,                                                                                           \
										complex_functions::VECOP<T>{},                                                                                                 \
										complex_functions::OP<T>{});                                                                                                   \
	}                                                                                                                                                                  \
                                                                                                                                                                       \
	template <class SignalT, class SignalU, std::enable_if_t<is_mutable_signal_v<SignalT> && is_same_domain_v<std::decay_t<SignalT>, std::decay_t<SignalU>>, int> = 0> \
	auto NAME(SignalT&& out, const SignalU& in, int, ...) {                                                                                                            \
                                                                                                                                                                       \
		return UnaryOperationVectorized(out.Data(), in.Data(), out.Length(), [](auto v) { return math_functions::FUNC(v); });                                          \
	}                                                                                                                                                                  \
                                                                                                                                                                       \
	template <class SignalT, class SignalU, std::enable_if_t<is_mutable_signal_v<SignalT> && is_same_domain_v<std::decay_t<SignalT>, std::decay_t<SignalU>>, int> = 0> \
	auto NAME(SignalT&& out, const SignalU& in) {                                                                                                                      \
		return NAME(std::forward<SignalT>(out), in, 0, typename signal_traits<SignalU>::type{});                                                                       \
	}

#define DSPBB_IMPL_COMPLEX_FUNCTION_1_PARAM(NAME, FUNC)                                          \
	template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0> \
	auto NAME(const SignalT& signal) {                                                           \
		using R = decltype(std::FUNC(std::declval<typename signal_traits<SignalT>::type>()));    \
		Signal<R, signal_traits<SignalT>::domain> r(signal.Size());                              \
		NAME(r, signal);                                                                         \
		return r;                                                                                \
	}

#define DSPBB_IMPL_COMPLEX_FUNCTION(NAME, VECOP, OP, FUNC)     \
	DSPBB_IMPL_COMPLEX_FUNCTION_2_PARAM(NAME, VECOP, OP, FUNC) \
	DSPBB_IMPL_COMPLEX_FUNCTION_1_PARAM(NAME, FUNC)

DSPBB_IMPL_COMPLEX_FUNCTION(Abs, AbsVec, Abs, abs)
DSPBB_IMPL_COMPLEX_FUNCTION(Arg, ArgVec, Arg, arg)
DSPBB_IMPL_COMPLEX_FUNCTION(Real, RealVec, Real, real)
DSPBB_IMPL_COMPLEX_FUNCTION(Imag, ImagVec, Imag, imag)

//------------------------------------------------------------------------------
// Exponential functions
//------------------------------------------------------------------------------

DSPBB_IMPL_FUNCTION(Log, log)
DSPBB_IMPL_FUNCTION(Log2, log2)
DSPBB_IMPL_FUNCTION(Log10, log10)
DSPBB_IMPL_FUNCTION(Exp, exp)

//------------------------------------------------------------------------------
// Polynomial functions
//------------------------------------------------------------------------------


template <class SignalT, class SignalU, std::enable_if_t<is_mutable_signal_v<SignalT> && is_same_domain_v<std::decay_t<SignalT>, std::decay_t<SignalU>>, int> = 0>
auto Pow(SignalT&& out, const SignalU& in, std::remove_const_t<typename std::decay_t<SignalU>::value_type> power) {
	return UnaryOperationVectorized(out.Data(), in.Data(), out.Length(), [power](auto v) { return math_functions::pow(v, power); });
}
template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto Pow(const SignalT& signal, std::remove_const_t<typename std::decay_t<SignalT>::value_type> power) {
	SignalT r(signal.Size());
	Pow(r, signal, power);
	return r;
}

DSPBB_IMPL_FUNCTION(Sqrt, sqrt)
DSPBB_IMPL_FUNCTION(Cbrt, cbrt)

//------------------------------------------------------------------------------
// Trigonometric functions
//------------------------------------------------------------------------------

DSPBB_IMPL_FUNCTION(Sin, sin)
DSPBB_IMPL_FUNCTION(Cos, cos)
DSPBB_IMPL_FUNCTION(Tan, tan)
DSPBB_IMPL_FUNCTION(Asin, asin)
DSPBB_IMPL_FUNCTION(Acos, acos)
DSPBB_IMPL_FUNCTION(Atan, atan)


//------------------------------------------------------------------------------
// Hyperbolic functions
//------------------------------------------------------------------------------

DSPBB_IMPL_FUNCTION(Sinh, sinh)
DSPBB_IMPL_FUNCTION(Cosh, cosh)
DSPBB_IMPL_FUNCTION(Tanh, tanh)
DSPBB_IMPL_FUNCTION(Asinh, asinh)
DSPBB_IMPL_FUNCTION(Acosh, acosh)
DSPBB_IMPL_FUNCTION(Atanh, atanh)


} // namespace dspbb