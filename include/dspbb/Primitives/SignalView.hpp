#pragma once

#include "../Utility/TypeTraits.hpp"
#include "Signal.hpp"


namespace dspbb {


template <class T, eSignalDomain Domain>
class BasicSignalView {
public:
	using iterator = T*;
	using const_iterator = const T*;
	using reverse_iterator = std::reverse_iterator<iterator>;
	using const_reverse_iterator = std::reverse_iterator<const_iterator>;
	using size_type = std::size_t;
	using value_type = T;

public:
	BasicSignalView() = default;
	BasicSignalView(BasicSignalView&&) noexcept = default;
	BasicSignalView(const BasicSignalView&) noexcept = default;
	BasicSignalView& operator=(BasicSignalView&&) noexcept = default;
	BasicSignalView& operator=(const BasicSignalView&) noexcept = default;

	explicit BasicSignalView(BasicSignal<std::remove_const_t<T>, Domain>& signal);

	template <class Q = T, std::enable_if_t<std::is_const_v<Q>, int> = 0>
	explicit BasicSignalView(const BasicSignal<std::remove_const_t<T>, Domain>& signal);

	template <class Q = T, std::enable_if_t<std::is_const_v<Q>, int> = 0>
	BasicSignalView(const BasicSignalView<std::remove_const_t<T>, Domain>& signal);

	template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter&>()), T&>, int> = 0>
	BasicSignalView(Iter first, Iter last);
	template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter&>()), T&>, int> = 0>
	BasicSignalView(Iter first, size_t size);

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

	BasicSignalView First(size_type n);
	BasicSignalView Last(size_type n);
	BasicSignalView SubSignal(size_type offset) const;
	BasicSignalView SubSignal(size_type offset, size_type count) const;

private:
	iterator first = nullptr, last = nullptr;
};

template <class T, eSignalDomain Domain>
BasicSignalView<T, Domain>::BasicSignalView(BasicSignal<std::remove_const_t<T>, Domain>& signal)
	: BasicSignalView(signal.begin(), signal.end()) {
}

template <class T, eSignalDomain Domain>
template <class Q, std::enable_if_t<std::is_const_v<Q>, int>>
BasicSignalView<T, Domain>::BasicSignalView(const BasicSignal<std::remove_const_t<T>, Domain>& signal)
	: BasicSignalView(signal.begin(), signal.end()) {
}

template <class T, eSignalDomain Domain>
template <class Q, std::enable_if_t<std::is_const_v<Q>, int>>
BasicSignalView<T, Domain>::BasicSignalView(const BasicSignalView<std::remove_const_t<T>, Domain>& signal)
	: BasicSignalView(signal.begin(), signal.end()) {
}

template <class T, eSignalDomain Domain>
template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter&>()), T&>, int>>
BasicSignalView<T, Domain>::BasicSignalView(Iter first, Iter last) {
	this->first = first != last ? std::addressof(*first) : nullptr;
	this->last = this->first + (last - first);
}

template <class T, eSignalDomain Domain>
template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter&>()), T&>, int>>
BasicSignalView<T, Domain>::BasicSignalView(Iter first, size_t size) {
	this->first = size != 0 ? std::addressof(*first) : nullptr;
	this->last = this->first + size;
}


template <class T, eSignalDomain Domain>
T& BasicSignalView<T, Domain>::Front() const { return *first; }

template <class T, eSignalDomain Domain>
T& BasicSignalView<T, Domain>::Back() const { return *(last - 1); }

template <class T, eSignalDomain Domain>
T& BasicSignalView<T, Domain>::operator[](size_type index) const { return first[index]; }

template <class T, eSignalDomain Domain>
T* BasicSignalView<T, Domain>::Data() const { return std::addressof(*first); }

template <class T, eSignalDomain Domain>
typename BasicSignalView<T, Domain>::size_type BasicSignalView<T, Domain>::Size() const { return size_type(last - first); }

template <class T, eSignalDomain Domain>
typename BasicSignalView<T, Domain>::size_type BasicSignalView<T, Domain>::Length() const { return Size(); }

template <class T, eSignalDomain Domain>
typename BasicSignalView<T, Domain>::size_type BasicSignalView<T, Domain>::SizeBytes() const { return sizeof(T) * Size(); }

template <class T, eSignalDomain Domain>
bool BasicSignalView<T, Domain>::Empty() const { return first == last; }

template <class T, eSignalDomain Domain>
BasicSignalView<T, Domain> BasicSignalView<T, Domain>::First(size_type n) {
	assert(n < Size());
	return { first, first + n };
}

template <class T, eSignalDomain Domain>
BasicSignalView<T, Domain> BasicSignalView<T, Domain>::Last(size_type n) {
	assert(n < Size());
	return { last - n, last };
}

template <class T, eSignalDomain Domain>
BasicSignalView<T, Domain> BasicSignalView<T, Domain>::SubSignal(size_type offset) const {
	assert(offset <= this->Size());
	return { this->first + offset, this->last };
}

template <class T, eSignalDomain Domain>
BasicSignalView<T, Domain> BasicSignalView<T, Domain>::SubSignal(size_type offset, size_type count) const {
	assert(offset <= this->Size());
	assert(offset + count <= this->Size());
	return { this->first + offset, this->first + offset + count };
}

// Helpers
template <class T, eSignalDomain Domain>
auto AsView(BasicSignal<T, Domain>& signal) -> BasicSignalView<T, Domain> {
	return BasicSignalView<T, Domain>{ signal };
}

template <class T, eSignalDomain Domain>
auto AsView(const BasicSignal<T, Domain>& signal) -> BasicSignalView<const T, Domain> {
	return BasicSignalView<const T, Domain>{ signal };
}

template <class T, eSignalDomain Domain>
auto AsView(BasicSignalView<T, Domain> view) -> BasicSignalView<T, Domain> {
	return view;
}

template <class T, eSignalDomain Domain>
auto AsView(BasicSignalView<const T, Domain> view) -> BasicSignalView<const T, Domain> {
	return view;
}

template <class T, eSignalDomain Domain>
auto AsConstView(const BasicSignal<T, Domain>& signal) -> BasicSignalView<const T, Domain> {
	return BasicSignalView<const T, Domain>{ signal };
}

template <class T, eSignalDomain Domain>
auto AsConstView(BasicSignalView<T, Domain> view) -> BasicSignalView<const T, Domain> {
	return view;
}

template <class T, eSignalDomain Domain>
auto AsConstView(BasicSignalView<const T, Domain> view) -> BasicSignalView<const T, Domain> {
	return view;
}

template <eSignalDomain Domain, class Iter>
auto AsView(Iter first, Iter last) {
	using T = typename std::iterator_traits<Iter>::value_type;
	return BasicSignalView<T, Domain>{ first, last };
}

template <eSignalDomain Domain, class Iter>
auto AsView(Iter first, size_t size) {
	using T = typename std::iterator_traits<Iter>::value_type;
	return BasicSignalView<T, Domain>{ first, size };
}

template <eSignalDomain Domain, class Iter>
auto AsConstView(Iter first, Iter last) {
	using T = typename std::iterator_traits<Iter>::value_type;
	return BasicSignalView<const T, Domain>{ first, last };
}

template <eSignalDomain Domain, class Iter>
auto AsConstView(Iter first, size_t size) {
	using T = typename std::iterator_traits<Iter>::value_type;
	return BasicSignalView<const T, Domain>{ first, size };
}


template <class T>
using SignalView = BasicSignalView<T, eSignalDomain::TIME>;
template <class T>
using SpectrumView = BasicSignalView<T, eSignalDomain::FREQUENCY>;
template <class T>
using CepstrumView = BasicSignalView<T, eSignalDomain::QUEFRENCY>;

using SignalViewF = SignalView<float>;
using SignalViewCF = SignalView<std::complex<float>>;

using SpectrumViewCF = BasicSignalView<std::complex<float>, eSignalDomain::FREQUENCY>;
using SpectrumViewF = BasicSignalView<float, eSignalDomain::FREQUENCY>;

using CepstrumViewCF = BasicSignalView<std::complex<float>, eSignalDomain::QUEFRENCY>;
using CepstrumViewF = BasicSignalView<float, eSignalDomain::QUEFRENCY>;


} // namespace dspbb

#include "SignalArithmetic.hpp"