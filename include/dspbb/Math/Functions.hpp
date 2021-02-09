#pragma once

#include "../Primitives/Signal.hpp"
#include "../Primitives/SignalView.hpp"
#include "../Utility/Algorithm.hpp"

#include <complex>
#include <type_traits>


namespace dspbb {


//------------------------------------------------------------------------------
// Complex number functions
//------------------------------------------------------------------------------


template <class SignalT>
auto Abs(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::abs(v); });
}

template <class SignalT, std::enable_if_t<is_complex_v<typename std::decay_t<SignalT>::value_type>, int> = 0>
auto Arg(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::arg(v); });
}

template <class SignalT, std::enable_if_t<!is_complex_v<typename std::decay_t<SignalT>::value_type>, int> = 0>
auto Real(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return v; });
}

template <class SignalT, std::enable_if_t<is_complex_v<typename std::decay_t<SignalT>::value_type>, int> = 0>
auto Real(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::real(v); });
}

template <class SignalT, std::enable_if_t<is_complex_v<typename std::decay_t<SignalT>::value_type>, int> = 0>
auto Imag(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::imag(v); });
}


//------------------------------------------------------------------------------
// Exponential functions
//------------------------------------------------------------------------------

template <class SignalT>
auto Log(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::log(v); });
}

template <class SignalT>
auto Log2(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::log2(v); });
}

template <class SignalT>
auto Log10(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::log10(v); });
}

template <class SignalT>
auto Exp(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::exp(v); });
}


//------------------------------------------------------------------------------
// Polynomial functions
//------------------------------------------------------------------------------


template <class SignalT>
auto Pow(SignalT&& signal, typename std::decay_t<SignalT>::value_type power) {
	return Apply(
		std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v, typename std::decay_t<SignalT>::value_type power) { return std::pow(v, power); }, power);
}

template <class SignalT>
auto Sqrt(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::sqrt(v); });
}

template <class SignalT>
auto Cbrt(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::cbrt(v); });
}


//------------------------------------------------------------------------------
// Trigonometric functions
//------------------------------------------------------------------------------

template <class SignalT>
auto Sin(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::sin(v); });
}

template <class SignalT>
auto Cos(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::cos(v); });
}


template <class SignalT>
auto Tan(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::tan(v); });
}

template <class SignalT>
auto Asin(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::asin(v); });
}

template <class SignalT>
auto Acos(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::acos(v); });
}

template <class SignalT>
auto Atan(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::atan(v); });
}


//------------------------------------------------------------------------------
// Hyperbolic functions
//------------------------------------------------------------------------------

template <class SignalT>
auto Sinh(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::sinh(v); });
}

template <class SignalT>
auto Cosh(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::cosh(v); });
}


template <class SignalT>
auto Tanh(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::tanh(v); });
}

template <class SignalT>
auto Asinh(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::asinh(v); });
}

template <class SignalT>
auto Acosh(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::acosh(v); });
}

template <class SignalT>
auto Atanh(SignalT&& signal) {
	return Apply(std::forward<SignalT>(signal), [](typename std::decay_t<SignalT>::value_type v) { return std::atanh(v); });
}


} // namespace dspbb