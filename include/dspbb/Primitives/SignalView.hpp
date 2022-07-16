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

	BasicSignalView(BasicSignal<std::remove_const_t<T>, Domain>& signal);

	template <class Q = T, std::enable_if_t<std::is_const_v<Q>, int> = 0>
	BasicSignalView(const BasicSignal<std::remove_const_t<T>, Domain>& signal);

	template <class Q = T, std::enable_if_t<std::is_const_v<Q>, int> = 0>
	BasicSignalView(const BasicSignalView<std::remove_const_t<T>, Domain>& signal);

	template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter&>()), T&>, int> = 0>
	BasicSignalView(Iter first, Iter last);
	template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter&>()), T&>, int> = 0>
	BasicSignalView(Iter first, size_t size);

	iterator begin() const { return m_first; }
	const_iterator cbegin() const { return m_first; }
	iterator end() const { return m_last; }
	const_iterator cend() const { return m_last; }
	reverse_iterator rbegin() const { return reverse_iterator{ m_last }; }
	const_reverse_iterator crbegin() const { return const_reverse_iterator{ m_last }; }
	reverse_iterator rend() const { return reverse_iterator{ m_first }; }
	const_reverse_iterator crend() const { return const_reverse_iterator{ m_first }; }

	T& front() const;
	T& back() const;
	T& operator[](size_type index) const;
	T* data() const;

	size_type size() const;
	size_type size_bytes() const;
	bool empty() const;

	BasicSignalView first(size_type n);
	BasicSignalView last(size_type n);
	BasicSignalView subsignal(size_type offset) const;
	BasicSignalView subsignal(size_type offset, size_type count) const;

private:
	iterator m_first = nullptr, m_last = nullptr;
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
	this->m_first = first != last ? std::addressof(*first) : nullptr;
	this->m_last = this->m_first + (last - first);
}

template <class T, eSignalDomain Domain>
template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter&>()), T&>, int>>
BasicSignalView<T, Domain>::BasicSignalView(Iter first, size_t size) {
	this->m_first = size != 0 ? std::addressof(*first) : nullptr;
	this->m_last = this->m_first + size;
}


template <class T, eSignalDomain Domain>
T& BasicSignalView<T, Domain>::front() const { return *m_first; }

template <class T, eSignalDomain Domain>
T& BasicSignalView<T, Domain>::back() const { return *(m_last - 1); }

template <class T, eSignalDomain Domain>
T& BasicSignalView<T, Domain>::operator[](size_type index) const { return m_first[index]; }

template <class T, eSignalDomain Domain>
T* BasicSignalView<T, Domain>::data() const { return std::addressof(*m_first); }

template <class T, eSignalDomain Domain>
typename BasicSignalView<T, Domain>::size_type BasicSignalView<T, Domain>::size() const { return size_type(m_last - m_first); }

template <class T, eSignalDomain Domain>
typename BasicSignalView<T, Domain>::size_type BasicSignalView<T, Domain>::size_bytes() const { return sizeof(T) * size(); }

template <class T, eSignalDomain Domain>
bool BasicSignalView<T, Domain>::empty() const { return m_first == m_last; }

template <class T, eSignalDomain Domain>
BasicSignalView<T, Domain> BasicSignalView<T, Domain>::first(size_type n) {
	assert(n <= size());
	return { m_first, m_first + n };
}

template <class T, eSignalDomain Domain>
BasicSignalView<T, Domain> BasicSignalView<T, Domain>::last(size_type n) {
	assert(n <= size());
	return { m_last - n, m_last };
}

template <class T, eSignalDomain Domain>
BasicSignalView<T, Domain> BasicSignalView<T, Domain>::subsignal(size_type offset) const {
	assert(offset <= this->size());
	return { this->m_first + offset, this->m_last };
}

template <class T, eSignalDomain Domain>
BasicSignalView<T, Domain> BasicSignalView<T, Domain>::subsignal(size_type offset, size_type count) const {
	assert(offset <= this->size());
	assert(offset + count <= this->size());
	return { this->m_first + offset, this->m_first + offset + count };
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