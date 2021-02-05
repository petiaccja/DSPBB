#pragma once

#include <complex>
#include <type_traits>
#include "TypeTraits.hpp"

namespace dspbb {


template <class T, bool Pred>
struct add_const_conditional {
	using type = std::conditional_t<Pred, std::add_const_t<T>, T>;
};

template <class T, bool Pred>
using add_const_conditional_t = typename add_const_conditional<T, Pred>::type;

template <class T, bool Pred>
struct add_volatile_conditional {
	using type = std::conditional_t<Pred, std::add_volatile_t<T>, T>;
};

template <class T, bool Pred>
using add_volatile_conditional_t = typename add_volatile_conditional<T, Pred>::type;


template <class T, bool Const, bool Volatile>
struct add_cv_conditional {
	using type = add_volatile_conditional_t<add_const_conditional_t<T, Const>, Volatile>;
};

template <class T, bool Const, bool Volatile>
using add_cv_conditional_t = typename add_cv_conditional<T, Const, Volatile>::type;


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