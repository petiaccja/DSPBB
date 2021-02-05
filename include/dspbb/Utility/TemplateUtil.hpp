#pragma once

#include <type_traits>

namespace enl {

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


}