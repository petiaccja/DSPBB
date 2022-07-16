#pragma once

#include "../Kernels/Functors.hpp"
#include "../Kernels/Numeric.hpp"
#include "SignalTraits.hpp"

#include <functional>


namespace dspbb {


//------------------------------------------------------------------------------
// Helpers.
//------------------------------------------------------------------------------

template <class SignalR, class SignalT, class SignalU>
void CheckSizes(const SignalR& r, const SignalT& a, const SignalU& b) {
	assert(r.size() == a.size());
	assert(r.size() == b.size());
	if (r.size() != a.size() || r.size() != b.size()) {
		throw std::invalid_argument("All input vectors must be the same size.");
	}
}

template <class SignalR, class SignalT>
void CheckSizes(const SignalR& r, const SignalT& a) {
	assert(r.size() == a.size());
	if (r.size() != a.size()) {
		throw std::invalid_argument("All input vectors must be the same size.");
	}
}

//------------------------------------------------------------------------------
// Three operand functions.
//------------------------------------------------------------------------------

template <class SignalR, class SignalT, class SignalU>
auto Multiply(SignalR&& r, const SignalT& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT, SignalU>, void> {
	CheckSizes(r, a, b);
	kernels::Transform(a.begin(), a.end(), b.begin(), r.begin(), std::multiplies{});
}

template <class SignalR, class T, class SignalU>
auto Multiply(SignalR&& r, const T& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalU> && !is_signal_like_v<T>, void> {
	CheckSizes(r, b);
	kernels::Transform(b.begin(), b.end(), r.begin(), multiplies_scalar_left{ a });
}

template <class SignalR, class SignalT, class U>
auto Multiply(SignalR&& r, const SignalT& a, const U& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT> && !is_signal_like_v<U>, void> {
	CheckSizes(r, a);
	kernels::Transform(a.begin(), a.end(), r.begin(), multiplies_scalar_right{ b });
}


template <class SignalR, class SignalT, class SignalU>
auto Divide(SignalR&& r, const SignalT& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT, SignalU>, void> {
	CheckSizes(r, a, b);
	kernels::Transform(a.begin(), a.end(), b.begin(), r.begin(), std::divides{});
}

template <class SignalR, class T, class SignalU>
auto Divide(SignalR&& r, const T& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalU> && !is_signal_like_v<T>, void> {
	CheckSizes(r, b);
	kernels::Transform(b.begin(), b.end(), r.begin(), divides_scalar_left{ a });
}

template <class SignalR, class SignalT, class U>
auto Divide(SignalR&& r, const SignalT& a, const U& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT> && !is_signal_like_v<U>, void> {
	CheckSizes(r, a);
	kernels::Transform(a.begin(), a.end(), r.begin(), divides_scalar_right{ b });
}


template <class SignalR, class SignalT, class SignalU>
auto Add(SignalR&& r, const SignalT& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT, SignalU>, void> {
	CheckSizes(r, a, b);
	kernels::Transform(a.begin(), a.end(), b.begin(), r.begin(), std::plus{});
}

template <class SignalR, class T, class SignalU>
auto Add(SignalR&& r, const T& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalU> && !is_signal_like_v<T>, void> {
	CheckSizes(r, b);
	kernels::Transform(b.begin(), b.end(), r.begin(), plus_scalar_left{ a });
}

template <class SignalR, class SignalT, class U>
auto Add(SignalR&& r, const SignalT& a, const U& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT> && !is_signal_like_v<U>, void> {
	CheckSizes(r, a);
	kernels::Transform(a.begin(), a.end(), r.begin(), plus_scalar_right{ b });
}


template <class SignalR, class SignalT, class SignalU>
auto Subtract(SignalR&& r, const SignalT& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT, SignalU>, void> {
	CheckSizes(r, a, b);
	kernels::Transform(a.begin(), a.end(), b.begin(), r.begin(), std::minus{});
}

template <class SignalR, class T, class SignalU>
auto Subtract(SignalR&& r, const T& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalU> && !is_signal_like_v<T>, void> {
	CheckSizes(r, b);
	kernels::Transform(b.begin(), b.end(), r.begin(), minus_scalar_left{ a });
}

template <class SignalR, class SignalT, class U>
auto Subtract(SignalR&& r, const SignalT& a, const U& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT> && !is_signal_like_v<U>, void> {
	CheckSizes(r, a);
	kernels::Transform(a.begin(), a.end(), r.begin(), minus_scalar_right{ b });
}



//------------------------------------------------------------------------------
// Overloaded operators.
//------------------------------------------------------------------------------

//--------------------------------------
// Vector-vector
//--------------------------------------

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto operator*(const SignalT& a, const SignalU& b) {
	using R = decltype(std::declval<typename signal_traits<SignalT>::type>() * std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalT>::domain;
	BasicSignal<R, Domain> r(a.size());
	Multiply<BasicSignal<R, Domain>&, SignalT, SignalU>(r, a, b);
	return r;
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto operator/(const SignalT& a, const SignalU& b) {
	using R = decltype(std::declval<typename signal_traits<SignalT>::type>() / std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalT>::domain;
	BasicSignal<R, Domain> r(a.size());
	Divide(r, a, b);
	return r;
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto operator+(const SignalT& a, const SignalU& b) {
	using R = decltype(std::declval<typename signal_traits<SignalT>::type>() + std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalT>::domain;
	BasicSignal<R, Domain> r(a.size());
	Add(r, a, b);
	return r;
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto operator-(const SignalT& a, const SignalU& b) {
	using R = decltype(std::declval<typename signal_traits<SignalT>::type>() - std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalT>::domain;
	BasicSignal<R, Domain> r(a.size());
	Subtract(r, a, b);
	return r;
}


//--------------------------------------
// Vector-scalar
//--------------------------------------

template <class SignalT, class U, std::enable_if_t<is_signal_like_v<SignalT> && !is_signal_like_v<U>, int> = 0>
auto operator*(const SignalT& a, const U& b) {
	using R = decltype(std::declval<typename signal_traits<SignalT>::type>() * std::declval<U>());
	constexpr auto Domain = signal_traits<SignalT>::domain;
	BasicSignal<R, Domain> r(a.size());
	Multiply(r, a, b);
	return r;
}

template <class SignalT, class U, std::enable_if_t<is_signal_like_v<SignalT> && !is_signal_like_v<U>, int> = 0>
auto operator/(const SignalT& a, const U& b) {
	using R = decltype(std::declval<typename signal_traits<SignalT>::type>() / std::declval<U>());
	constexpr auto Domain = signal_traits<SignalT>::domain;
	BasicSignal<R, Domain> r(a.size());
	Divide(r, a, b);
	return r;
}

template <class SignalT, class U, std::enable_if_t<is_signal_like_v<SignalT> && !is_signal_like_v<U>, int> = 0>
auto operator+(const SignalT& a, const U& b) {
	using R = decltype(std::declval<typename signal_traits<SignalT>::type>() + std::declval<U>());
	constexpr auto Domain = signal_traits<SignalT>::domain;
	BasicSignal<R, Domain> r(a.size());
	Add(r, a, b);
	return r;
}

template <class SignalT, class U, std::enable_if_t<is_signal_like_v<SignalT> && !is_signal_like_v<U>, int> = 0>
auto operator-(const SignalT& a, const U& b) {
	using R = decltype(std::declval<typename signal_traits<SignalT>::type>() - std::declval<U>());
	constexpr auto Domain = signal_traits<SignalT>::domain;
	BasicSignal<R, Domain> r(a.size());
	Subtract(r, a, b);
	return r;
}


template <class T, class SignalU, std::enable_if_t<!is_signal_like_v<T> && is_signal_like_v<std::decay_t<SignalU>>, int> = 0>
auto operator*(const T& a, const SignalU& b) {
	using R = decltype(std::declval<T>() * std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalU>::domain;
	BasicSignal<R, Domain> r(b.size());
	Multiply(r, a, b);
	return r;
}

template <class T, class SignalU, std::enable_if_t<!is_signal_like_v<T> && is_signal_like_v<std::decay_t<SignalU>>, int> = 0>
auto operator/(const T& a, const SignalU& b) {
	using R = decltype(std::declval<T>() / std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalU>::domain;
	BasicSignal<R, Domain> r(b.size());
	Divide(r, a, b);
	return r;
}

template <class T, class SignalU, std::enable_if_t<!is_signal_like_v<T> && is_signal_like_v<std::decay_t<SignalU>>, int> = 0>
auto operator+(const T& a, const SignalU& b) {
	using R = decltype(std::declval<T>() + std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalU>::domain;
	BasicSignal<R, Domain> r(b.size());
	Add(r, a, b);
	return r;
}

template <class T, class SignalU, std::enable_if_t<!is_signal_like_v<T> && is_signal_like_v<std::decay_t<SignalU>>, int> = 0>
auto operator-(const T& a, const SignalU& b) {
	using R = decltype(std::declval<T>() - std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalU>::domain;
	BasicSignal<R, Domain> r(b.size());
	Subtract(r, a, b);
	return r;
}


//--------------------------------------
// Compound vector-vector
//--------------------------------------

template <class SignalT, class SignalU>
auto operator*=(SignalT&& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalT&> && is_same_domain_v<SignalT, SignalU>, SignalT&> {
	Multiply(a, a, b);
	return a;
}

template <class SignalT, class SignalU>
auto operator/=(SignalT&& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalT&> && is_same_domain_v<SignalT, SignalU>, SignalT&> {
	Divide(a, a, b);
	return a;
}

template <class SignalT, class SignalU>
auto operator+=(SignalT&& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalT&> && is_same_domain_v<SignalT, SignalU>, SignalT&> {
	Add(a, a, b);
	return a;
}

template <class SignalT, class SignalU>
auto operator-=(SignalT&& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalT&> && is_same_domain_v<SignalT, SignalU>, SignalT&> {
	Subtract(a, a, b);
	return a;
}


//--------------------------------------
// Compound vector-scalar
//--------------------------------------

template <class SignalT, class U>
auto operator*=(SignalT&& a, const U& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalT&> && is_signal_like_v<std::decay_t<SignalT>> && !is_signal_like_v<U>, SignalT&> {
	Multiply(a, a, b);
	return a;
}

template <class SignalT, class U>
auto operator/=(SignalT&& a, const U& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalT&> && is_signal_like_v<std::decay_t<SignalT>> && !is_signal_like_v<U>, SignalT&> {
	Divide(a, a, b);
	return a;
}

template <class SignalT, class U>
auto operator+=(SignalT&& a, const U& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalT&> && is_signal_like_v<std::decay_t<SignalT>> && !is_signal_like_v<U>, SignalT&> {
	Add(a, a, b);
	return a;
}

template <class SignalT, class U>
auto operator-=(SignalT&& a, const U& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalT&> && is_signal_like_v<std::decay_t<SignalT>> && !is_signal_like_v<U>, SignalT&> {
	Subtract(a, a, b);
	return a;
}


} // namespace dspbb