#pragma once

#include "../Utility/TypeTraits.hpp"
#include "Signal.hpp"


namespace dspbb {


template <class T, eSignalDomain Domain>
class SignalView {
public:
	using iterator = T*;
	using const_iterator = const T*;
	using reverse_iterator = std::reverse_iterator<iterator>;
	using const_reverse_iterator = std::reverse_iterator<const_iterator>;
	using size_type = std::size_t;
	using value_type = T;

public:
	SignalView() = default;
	SignalView(SignalView&&) noexcept = default;
	SignalView(const SignalView&) noexcept = default;
	SignalView& operator=(SignalView&&) noexcept = default;
	SignalView& operator=(const SignalView&) noexcept = default;

	explicit SignalView(BasicSignal<std::remove_const_t<T>, Domain>& signal);

	template <class Q = T, std::enable_if_t<std::is_const_v<Q>, int> = 0>
	explicit SignalView(const BasicSignal<std::remove_const_t<T>, Domain>& signal);

	template <class Q = T, std::enable_if_t<std::is_const_v<Q>, int> = 0>
	SignalView(const SignalView<std::remove_const_t<T>, Domain>& signal);

	template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter&>()), T&>, int> = 0>
	SignalView(Iter first, Iter last);
	template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter&>()), T&>, int> = 0>
	SignalView(Iter first, size_t size);

	iterator begin() const { return first; }
	const_iterator cbegin() const { return first; }
	iterator end() const { return last; }
	const_iterator cend() const { return last; }
	reverse_iterator rbegin() const { return reverse_iterator{ last }; }
	const_reverse_iterator crbegin() const { return const_reverse_iterator{ last }; }
	reverse_iterator rend() const { return reverse_iterator{ first }; }
	const_reverse_iterator crend() const { return const_reverse_iterator{ first }; }

	T& Front() const;
	T& Back() const;
	T& operator[](size_type index) const;
	T* Data() const;

	size_type Size() const;
	size_type Length() const;
	size_type SizeBytes() const;
	bool Empty() const;

	SignalView First(size_type n);
	SignalView Last(size_type n);
	SignalView SubSignal(size_type offset) const;
	SignalView SubSignal(size_type offset, size_type count) const;

private:
	iterator first = nullptr, last = nullptr;
};

template <class T, eSignalDomain Domain>
SignalView<T, Domain>::SignalView(BasicSignal<std::remove_const_t<T>, Domain>& signal)
	: SignalView(signal.begin(), signal.end()) {
}

template <class T, eSignalDomain Domain>
template <class Q, std::enable_if_t<std::is_const_v<Q>, int>>
SignalView<T, Domain>::SignalView(const BasicSignal<std::remove_const_t<T>, Domain>& signal)
	: SignalView(signal.begin(), signal.end()) {
}

template <class T, eSignalDomain Domain>
template <class Q, std::enable_if_t<std::is_const_v<Q>, int>>
SignalView<T, Domain>::SignalView(const SignalView<std::remove_const_t<T>, Domain>& signal)
	: SignalView(signal.begin(), signal.end()) {
}

template <class T, eSignalDomain Domain>
template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter&>()), T&>, int>>
SignalView<T, Domain>::SignalView(Iter first, Iter last) {
	this->first = first != last ? std::addressof(*first) : nullptr;
	this->last = this->first + (last - first);
}

template <class T, eSignalDomain Domain>
template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter&>()), T&>, int>>
SignalView<T, Domain>::SignalView(Iter first, size_t size) {
	this->first = size != 0 ? std::addressof(*first) : nullptr;
	this->last = this->first + size;
}


template <class T, eSignalDomain Domain>
T& SignalView<T, Domain>::Front() const { return *first; }

template <class T, eSignalDomain Domain>
T& SignalView<T, Domain>::Back() const { return *(last - 1); }

template <class T, eSignalDomain Domain>
T& SignalView<T, Domain>::operator[](size_type index) const { return first[index]; }

template <class T, eSignalDomain Domain>
T* SignalView<T, Domain>::Data() const { return std::addressof(*first); }

template <class T, eSignalDomain Domain>
typename SignalView<T, Domain>::size_type SignalView<T, Domain>::Size() const { return size_type(last - first); }

template <class T, eSignalDomain Domain>
typename SignalView<T, Domain>::size_type SignalView<T, Domain>::Length() const { return Size(); }

template <class T, eSignalDomain Domain>
typename SignalView<T, Domain>::size_type SignalView<T, Domain>::SizeBytes() const { return sizeof(T) * Size(); }

template <class T, eSignalDomain Domain>
bool SignalView<T, Domain>::Empty() const { return first == last; }

template <class T, eSignalDomain Domain>
SignalView<T, Domain> SignalView<T, Domain>::First(size_type n) {
	assert(n < Size());
	return { first, first + n };
}

template <class T, eSignalDomain Domain>
SignalView<T, Domain> SignalView<T, Domain>::Last(size_type n) {
	assert(n < Size());
	return { last - n, last };
}

template <class T, eSignalDomain Domain>
SignalView<T, Domain> SignalView<T, Domain>::SubSignal(size_type offset) const {
	assert(offset <= this->Size());
	return { this->first + offset, this->last };
}

template <class T, eSignalDomain Domain>
SignalView<T, Domain> SignalView<T, Domain>::SubSignal(size_type offset, size_type count) const {
	assert(offset <= this->Size());
	assert(offset + count <= this->Size());
	return { this->first + offset, this->first + offset + count };
}

// Helpers
template <class T, eSignalDomain Domain>
auto AsView(BasicSignal<T, Domain>& signal) -> SignalView<T, Domain> {
	return SignalView<T, Domain>{ signal };
}

template <class T, eSignalDomain Domain>
auto AsView(const BasicSignal<T, Domain>& signal) -> SignalView<const T, Domain> {
	return SignalView<const T, Domain>{ signal };
}

template <class T, eSignalDomain Domain>
auto AsView(SignalView<T, Domain> view) -> SignalView<T, Domain> {
	return view;
}

template <class T, eSignalDomain Domain>
auto AsView(SignalView<const T, Domain> view) -> SignalView<const T, Domain> {
	return view;
}

template <class T, eSignalDomain Domain>
auto AsConstView(const BasicSignal<T, Domain>& signal) -> SignalView<const T, Domain> {
	return SignalView<const T, Domain>{ signal };
}

template <class T, eSignalDomain Domain>
auto AsConstView(SignalView<T, Domain> view) -> SignalView<const T, Domain> {
	return view;
}

template <class T, eSignalDomain Domain>
auto AsConstView(SignalView<const T, Domain> view) -> SignalView<const T, Domain> {
	return view;
}

template <eSignalDomain Domain, class Iter>
auto AsView(Iter first, Iter last) {
	using T = typename std::iterator_traits<Iter>::value_type;
	return SignalView<T, Domain>{ first, last };
}

template <eSignalDomain Domain, class Iter>
auto AsView(Iter first, size_t size) {
	using T = typename std::iterator_traits<Iter>::value_type;
	return SignalView<T, Domain>{ first, size };
}

template <eSignalDomain Domain, class Iter>
auto AsConstView(Iter first, Iter last) {
	using T = typename std::iterator_traits<Iter>::value_type;
	return SignalView<const T, Domain>{ first, last };
}

template <eSignalDomain Domain, class Iter>
auto AsConstView(Iter first, size_t size) {
	using T = typename std::iterator_traits<Iter>::value_type;
	return SignalView<const T, Domain>{ first, size };
}


template <class T>
using TimeSignalView = SignalView<T, eSignalDomain::TIME>;
template <class T>
using SpectrumView = SignalView<T, eSignalDomain::FREQUENCY>;
template <class T>
using CepstrumView = SignalView<T, eSignalDomain::QUEFRENCY>;

using TimeSignalViewF = TimeSignalView<float>;
using TimeSignalViewCF = TimeSignalView<std::complex<float>>;

using SpectrumViewCF = SignalView<std::complex<float>, eSignalDomain::FREQUENCY>;
using SpectrumViewF = SignalView<float, eSignalDomain::FREQUENCY>;

using CepstrumViewCF = SignalView<std::complex<float>, eSignalDomain::QUEFRENCY>;
using CepstrumViewF = SignalView<float, eSignalDomain::QUEFRENCY>;


} // namespace dspbb

#include "SignalArithmetic.hpp"