#pragma once

#include "TypeTraits.hpp"

#include <complex>
#include <iterator>
#include <type_traits>

#pragma warning(push)
#pragma warning(disable : 4244)

namespace dspbb {

template <class T>
struct type_identity_cpp17 {
	using type = T;
};


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
	using type = add_cv_conditional_t<typename impl::remove_complex_h<std::remove_cv_t<T>>::type, std::is_const_v<T>, std::is_volatile_v<T>>;
};

template <class T>
using remove_complex_t = typename remove_complex<T>::type;


namespace impl {
	template <class VoidT, class T, class U>
	struct multiplies_result_helper {};

	template <class T, class U>
	struct multiplies_result_helper<std::void_t<decltype(std::declval<T>() * std::declval<U>())>, T, U> {
		using type = decltype(std::declval<T>() * std::declval<U>());
	};
} // namespace impl


template <class T, class U>
struct multiplies_result : impl::multiplies_result_helper<void, T, U> {};

template <class T, class U>
using multiplies_result_t = typename multiplies_result<T, U>::type;


namespace impl {
	template <class VoidT, class T, class U>
	struct plus_result_helper {};

	template <class T, class U>
	struct plus_result_helper<std::void_t<decltype(std::declval<T>() + std::declval<U>())>, T, U> {
		using type = decltype(std::declval<T>() + std::declval<U>());
	};
} // namespace impl

template <class T, class U>
struct plus_result : impl::plus_result_helper<void, T, U> {};

template <class T, class U>
using plus_result_t = typename plus_result<T, U>::type;


namespace impl {
	template <class VoidT, class T, class U>
	struct divides_result_helper {};

	template <class T, class U>
	struct divides_result_helper<std::void_t<decltype(std::declval<T>() / std::declval<U>())>, T, U> {
		using type = decltype(std::declval<T>() / std::declval<U>());
	};
} // namespace impl

template <class T, class U>
struct divides_result : impl::divides_result_helper<void, T, U> {};

template <class T, class U>
using divides_result_t = typename divides_result<T, U>::type;


namespace impl {
	template <class VoidT, class T, class U>
	struct minus_result_helper {};

	template <class T, class U>
	struct minus_result_helper<std::void_t<decltype(std::declval<T>() - std::declval<U>())>, T, U> {
		using type = decltype(std::declval<T>() - std::declval<U>());
	};
} // namespace impl

template <class T, class U>
struct minus_result : impl::minus_result_helper<void, T, U> {};

template <class T, class U>
using minus_result_t = typename minus_result<T, U>::type;


template <class Iter>
struct is_random_access_iterator : std::is_same<std::random_access_iterator_tag, typename std::iterator_traits<Iter>::iterator_category> {};

template <class Iter>
constexpr bool is_random_access_iterator_v = is_random_access_iterator<Iter>::value;


} // namespace dspbb

#pragma warning(pop)