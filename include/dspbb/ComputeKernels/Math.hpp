#pragma once

#include <cmath>
#ifdef _MSC_VER
	#pragma warning(push)
	#pragma warning(disable : 4800 4244)
#endif
#include <xsimd/xsimd.hpp>
#ifdef _MSC_VER
	#pragma warning(pop)
#endif

namespace dspbb::kernels {
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
	using std::imag;
	using std::real;

	using xsimd::abs;
	using xsimd::arg;
	using xsimd::imag;
	using xsimd::real;

	using std::conj;
	using xsimd::conj;

	// Erf & gamma
	using std::erf;
	using std::erfc;
	using std::lgamma;
	using std::tgamma;

	using xsimd::erf;
	using xsimd::erfc;
	using xsimd::lgamma;
	using xsimd::tgamma;

	// Misc
	namespace impl {
		template <class T, std::enable_if_t<(xsimd::simd_batch_traits<T>::size > 1), int> = 0>
		decltype(auto) min(const T& a, const T& b, int) {
			return xsimd::min(a, b);
		}
		template <class T, std::enable_if_t<(xsimd::simd_batch_traits<T>::size > 1), int> = 0>
		decltype(auto) max(const T& a, const T& b, int) {
			return xsimd::max(a, b);
		}
		template <class T>
		decltype(auto) min(const T& a, const T& b, ...) {
			return std::min(a, b);
		}
		template <class T>
		decltype(auto) max(const T& a, const T& b, ...) {
			return std::max(a, b);
		}
	} // namespace impl

	template <class T>
	decltype(auto) min(const T& a, const T& b) {
		return impl::min(a, b, 0);
	}
	template <class T>
	decltype(auto) max(const T& a, const T& b) {
		return impl::max(a, b, 0);
	}



} // namespace math_functions
} // namespace dspbb::kernels