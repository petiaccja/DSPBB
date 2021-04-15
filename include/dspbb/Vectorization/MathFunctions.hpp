#pragma once

#include <cmath>
#include <xsimd/xsimd.hpp>

namespace dspbb {
namespace math_functions {
	// Exponential
	using std::exp;
	using std::log;
	using std::log10;
	using std::log2;

	using xsimd::exp;
	using xsimd::log;
	using xsimd::log10;
	using xsimd::log2;

	// Polynomial
	using std::cbrt;
	using std::pow;
	using std::sqrt;

	using xsimd::cbrt;
	using xsimd::pow;
	using xsimd::sqrt;

	// Trigonometric
	using std::acos;
	using std::asin;
	using std::atan;
	using std::cos;
	using std::sin;
	using std::tan;

	using xsimd::acos;
	using xsimd::asin;
	using xsimd::atan;
	using xsimd::cos;
	using xsimd::sin;
	using xsimd::tan;

	// Hyperbolic
	using xsimd::acosh;
	using xsimd::asinh;
	using xsimd::atanh;
	using xsimd::cosh;
	using xsimd::sinh;
	using xsimd::tanh;

	using std::acosh;
	using std::asinh;
	using std::atanh;
	using std::cosh;
	using std::sinh;
	using std::tanh;

	// Complex
	using std::abs;
	using std::arg;
	using std::real;
	using std::imag;

	using xsimd::abs;
	using xsimd::arg;
	using xsimd::real;
	using xsimd::imag;

} // namespace math_functions
} // namespace dspbb