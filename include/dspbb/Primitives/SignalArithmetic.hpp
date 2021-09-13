#pragma once

#include "../Math/Arithmetic.hpp"
#include "SignalTraits.hpp"


namespace dspbb {


//------------------------------------------------------------------------------
// Helpers.
//------------------------------------------------------------------------------

template <class SignalR, class SignalT, class SignalU>
void CheckSizes(const SignalR& r, const SignalT& a, const SignalU& b) {
	assert(r.Size() == a.Size());
	assert(r.Size() == b.Size());
	if (r.Size() != a.Size() || r.Size() != b.Size()) {
		throw std::invalid_argument("All input vectors must be the same size.");
	}
}

template <class SignalR, class SignalT>
void CheckSizes(const SignalR& r, const SignalT& a) {
	assert(r.Size() == a.Size());
	if (r.Size() != a.Size()) {
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
	Multiply(r.Data(), a.Data(), b.Data(), r.Length());
}

template <class SignalR, class T, class SignalU>
auto Multiply(SignalR&& r, const T& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalU> && !is_signal_like_v<T>, void> {
	CheckSizes(r, b);
	Multiply(r.Data(), a, b.Data(), r.Length());
}

template <class SignalR, class SignalT, class U>
auto Multiply(SignalR&& r, const SignalT& a, const U& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT> && !is_signal_like_v<U>, void> {
	CheckSizes(r, a);
	Multiply(r.Data(), a.Data(), b, r.Length());
}


template <class SignalR, class SignalT, class SignalU>
auto Divide(SignalR&& r, const SignalT& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT, SignalU>, void> {
	CheckSizes(r, a, b);
	Divide(r.Data(), a.Data(), b.Data(), r.Length());
}

template <class SignalR, class T, class SignalU>
auto Divide(SignalR&& r, const T& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalU> && !is_signal_like_v<T>, void> {
	CheckSizes(r, b);
	Divide(r.Data(), a, b.Data(), r.Length());
}

template <class SignalR, class SignalT, class U>
auto Divide(SignalR&& r, const SignalT& a, const U& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT> && !is_signal_like_v<U>, void> {
	CheckSizes(r, a);
	Divide(r.Data(), a.Data(), b, r.Length());
}


template <class SignalR, class SignalT, class SignalU>
auto Add(SignalR&& r, const SignalT& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT, SignalU>, void> {
	CheckSizes(r, a, b);
	Add(r.Data(), a.Data(), b.Data(), r.Length());
}

template <class SignalR, class T, class SignalU>
auto Add(SignalR&& r, const T& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalU> && !is_signal_like_v<T>, void> {
	CheckSizes(r, b);
	Add(r.Data(), a, b.Data(), r.Length());
}

template <class SignalR, class SignalT, class U>
auto Add(SignalR&& r, const SignalT& a, const U& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT> && !is_signal_like_v<U>, void> {
	CheckSizes(r, a);
	Add(r.Data(), a.Data(), b, r.Length());
}


template <class SignalR, class SignalT, class SignalU>
auto Subtract(SignalR&& r, const SignalT& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT, SignalU>, void> {
	CheckSizes(r, a, b);
	Subtract(r.Data(), a.Data(), b.Data(), r.Length());
}

template <class SignalR, class T, class SignalU>
auto Subtract(SignalR&& r, const T& a, const SignalU& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalU> && !is_signal_like_v<T>, void> {
	CheckSizes(r, b);
	Subtract(r.Data(), a, b.Data(), r.Length());
}

template <class SignalR, class SignalT, class U>
auto Subtract(SignalR&& r, const SignalT& a, const U& b)
	-> std::enable_if_t<is_mutable_signal_v<SignalR&> && is_same_domain_v<SignalR, SignalT> && !is_signal_like_v<U>, void> {
	CheckSizes(r, a);
	Subtract(r.Data(), a.Data(), b, r.Length());
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
	Signal<R, Domain> r(a.Size());
	Multiply<Signal<R, Domain>&, SignalT, SignalU>(r, a, b);
	return r;
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto operator/(const SignalT& a, const SignalU& b) {
	using R = decltype(std::declval<typename signal_traits<SignalT>::type>() / std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalT>::domain;
	Signal<R, Domain> r(a.Size());
	Divide(r, a, b);
	return r;
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto operator+(const SignalT& a, const SignalU& b) {
	using R = decltype(std::declval<typename signal_traits<SignalT>::type>() + std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalT>::domain;
	Signal<R, Domain> r(a.Size());
	Add(r, a, b);
	return r;
}

template <class SignalT, class SignalU, std::enable_if_t<is_same_domain_v<SignalT, SignalU>, int> = 0>
auto operator-(const SignalT& a, const SignalU& b) {
	using R = decltype(std::declval<typename signal_traits<SignalT>::type>() - std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalT>::domain;
	Signal<R, Domain> r(a.Size());
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
	Signal<R, Domain> r(a.Size());
	Multiply(r, a, b);
	return r;
}

template <class SignalT, class U, std::enable_if_t<is_signal_like_v<SignalT> && !is_signal_like_v<U>, int> = 0>
auto operator/(const SignalT& a, const U& b) {
	using R = decltype(std::declval<typename signal_traits<SignalT>::type>() / std::declval<U>());
	constexpr auto Domain = signal_traits<SignalT>::domain;
	Signal<R, Domain> r(a.Size());
	Divide(r, a, b);
	return r;
}

template <class SignalT, class U, std::enable_if_t<is_signal_like_v<SignalT> && !is_signal_like_v<U>, int> = 0>
auto operator+(const SignalT& a, const U& b) {
	using R = decltype(std::declval<typename signal_traits<SignalT>::type>() + std::declval<U>());
	constexpr auto Domain = signal_traits<SignalT>::domain;
	Signal<R, Domain> r(a.Size());
	Add(r, a, b);
	return r;
}

template <class SignalT, class U, std::enable_if_t<is_signal_like_v<SignalT> && !is_signal_like_v<U>, int> = 0>
auto operator-(const SignalT& a, const U& b) {
	using R = decltype(std::declval<typename signal_traits<SignalT>::type>() - std::declval<U>());
	constexpr auto Domain = signal_traits<SignalT>::domain;
	Signal<R, Domain> r(a.Size());
	Subtract(r, a, b);
	return r;
}


template <class T, class SignalU, std::enable_if_t<!is_signal_like_v<T> && is_signal_like_v<std::decay_t<SignalU>>, int> = 0>
auto operator*(const T& a, const SignalU& b) {
	using R = decltype(std::declval<T>() * std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalU>::domain;
	Signal<R, Domain> r(b.Size());
	Multiply(r, a, b);
	return r;
}

template <class T, class SignalU, std::enable_if_t<!is_signal_like_v<T> && is_signal_like_v<std::decay_t<SignalU>>, int> = 0>
auto operator/(const T& a, const SignalU& b) {
	using R = decltype(std::declval<T>() / std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalU>::domain;
	Signal<R, Domain> r(b.Size());
	Divide(r, a, b);
	return r;
}

template <class T, class SignalU, std::enable_if_t<!is_signal_like_v<T> && is_signal_like_v<std::decay_t<SignalU>>, int> = 0>
auto operator+(const T& a, const SignalU& b) {
	using R = decltype(std::declval<T>() + std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalU>::domain;
	Signal<R, Domain> r(b.Size());
	Add(r, a, b);
	return r;
}

template <class T, class SignalU, std::enable_if_t<!is_signal_like_v<T> && is_signal_like_v<std::decay_t<SignalU>>, int> = 0>
auto operator-(const T& a, const SignalU& b) {
	using R = decltype(std::declval<T>() - std::declval<typename signal_traits<SignalU>::type>());
	constexpr auto Domain = signal_traits<SignalU>::domain;
	Signal<R, Domain> r(b.Size());
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