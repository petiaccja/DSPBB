#pragma once

#include "../ComputeKernels/VectorizedAlgorithms.hpp"
#include "../ComputeKernels/VectorizedMathFunctions.hpp"
#include "../Primitives/SignalTraits.hpp"

#include <complex>
#include <type_traits>


namespace dspbb {


#define DSPBB_IMPL_FUNCTION_2_PARAM(NAME, FUNC)                                                                                                                        \
	template <class SignalT, class SignalU, std::enable_if_t<is_mutable_signal_v<SignalT> && is_same_domain_v<std::decay_t<SignalT>, std::decay_t<SignalU>>, int> = 0> \
	auto NAME(SignalT&& out, const SignalU& in) {                                                                                                                      \
		return kernels::UnaryOperationVectorized(out.Data(), in.Data(), out.Length(), [](const auto& v) { return kernels::math_functions::FUNC(v); });                 \
	}

#define DSPBB_IMPL_FUNCTION_1_PARAM(NAME, FUNC)                                                  \
	template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0> \
	auto NAME(const SignalT& signal) {                                                           \
		using R = decltype(std::FUNC(std::declval<typename signal_traits<SignalT>::type>()));    \
		constexpr auto domain = signal_traits<SignalT>::domain;                                  \
		BasicSignal<R, domain> r(signal.Size());                                                      \
		NAME(r, signal);                                                                         \
		return r;                                                                                \
	}

#define DSPBB_IMPL_FUNCTION(NAME, FUNC)     \
	DSPBB_IMPL_FUNCTION_2_PARAM(NAME, FUNC) \
	DSPBB_IMPL_FUNCTION_1_PARAM(NAME, FUNC)


//------------------------------------------------------------------------------
// Complex number functions
//------------------------------------------------------------------------------

DSPBB_IMPL_FUNCTION(Abs, abs)
DSPBB_IMPL_FUNCTION(Arg, arg)
DSPBB_IMPL_FUNCTION(Real, real)
DSPBB_IMPL_FUNCTION(Imag, imag)
DSPBB_IMPL_FUNCTION(Conj, conj)

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
auto Pow(SignalT&& out, const SignalU& in, typename signal_traits<std::decay_t<SignalU>>::type power) {
	return kernels::UnaryOperationVectorized(out.Data(), in.Data(), out.Length(), [power](const auto& v) { return kernels::math_functions::pow(v, static_cast<std::decay_t<decltype(v)>>(power)); });
}
template <class SignalT, std::enable_if_t<is_signal_like_v<std::decay_t<SignalT>>, int> = 0>
auto Pow(const SignalT& signal, typename signal_traits<std::decay_t<SignalT>>::type power) {
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


//------------------------------------------------------------------------------
// Erf & gamma
//------------------------------------------------------------------------------

DSPBB_IMPL_FUNCTION(Erf, erf)
DSPBB_IMPL_FUNCTION(Erfc, erfc)
DSPBB_IMPL_FUNCTION(TGamma, tgamma)
DSPBB_IMPL_FUNCTION(LGamma, lgamma)

} // namespace dspbb