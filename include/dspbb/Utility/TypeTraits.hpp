#pragma once

#include <complex>
#include <type_traits>
#include "TemplateUtil.hpp"

namespace dspbb {

namespace impl {
	template <class T>
	struct is_complex_h : std::false_type {};

	template <class T>
	struct is_complex_h<std::complex<T>> : std::true_type {};

	template <class T>
	struct remove_complex_h {
		using type = std::decay_t<T>;
	};

	template <class T>
	struct remove_complex_h<std::complex<T>> {
		using type = T;
	};
} // namespace impl

template <class T>
struct is_complex : impl::is_complex_h<std::remove_cv_t<T>> {};

template <class T>
constexpr bool is_complex_v = is_complex<T>::value;


template <class T>
struct remove_complex : impl::remove_complex_h<std::remove_cv_t<T>> {
	using type = add_cv_conditional_t<typename impl::remove_complex_h<std::remove_cv_t<T>>::type, std::is_const<T>::value, std::is_volatile<T>::value>;
};

template <class T>
using remove_complex_t = typename remove_complex<T>::type;


} // namespace dspbb