#pragma once

#include "Definitions.hpp"

#include <iterator>
#include <span>


namespace dspbb {


template <class T, eSignalDomain Domain>
class BasicSignalView {
	using container = std::span<T>;

public:
	static constexpr bool is_const = std::is_const_v<T>;
	using element_type = typename container::element_type;
	using value_type = typename container::value_type;
	using pointer = typename container::pointer;
	using const_pointer = typename container::const_pointer;
	using reference = typename container::reference;
	using const_reference = typename container::const_reference;
	using size_type = typename container::size_type;

	using iterator = typename container::iterator;
	using reverse_iterator = typename container::reverse_iterator;

public:
	BasicSignalView() = default;
	BasicSignalView(BasicSignalView&&) noexcept = default;
	BasicSignalView(const BasicSignalView&) noexcept = default;
	BasicSignalView& operator=(BasicSignalView&&) noexcept = default;
	BasicSignalView& operator=(const BasicSignalView&) noexcept = default;

	explicit BasicSignalView(container c) : m_container(std::move(c)) {}

	BasicSignalView(BasicSignal<value_type, Domain>& signal)
		: m_container(signal.begin(), signal.end()) {}

	template <std::same_as<value_type> U>
	BasicSignalView(const BasicSignal<U, Domain>& signal) requires is_const
		: m_container(signal.begin(), signal.end()) {}

	template <std::same_as<value_type> U>
	BasicSignalView(const BasicSignalView<U, Domain>& signal) requires is_const
		: m_container(signal.begin(), signal.end()) {}

	template <std::contiguous_iterator Iter, std::sized_sentinel_for<Iter> End>
	BasicSignalView(Iter first, End last) : m_container(first, last) {}

	template <std::contiguous_iterator Iter>
	BasicSignalView(Iter first, size_t size) : m_container(first, size) {}

	reference front() const;
	reference back() const;
	reference operator[](size_type index) const;
	pointer data() const;

	size_type size() const;
	size_type size_bytes() const;
	bool empty() const;

	BasicSignalView first(size_type n);
	BasicSignalView last(size_type n);
	BasicSignalView subsignal(size_type offset) const;
	BasicSignalView subsignal(size_type offset, size_type count) const;

	iterator begin() const;
	iterator end() const;
	reverse_iterator rbegin() const;
	reverse_iterator rend() const;

private:
	container m_container;
};


template <class T, eSignalDomain Domain>
auto BasicSignalView<T, Domain>::front() const -> reference {
	return m_container.back();
}

template <class T, eSignalDomain Domain>
auto BasicSignalView<T, Domain>::back() const -> reference {
	return m_container.back();
}

template <class T, eSignalDomain Domain>
auto BasicSignalView<T, Domain>::operator[](size_type index) const -> reference {
	return m_container[index];
}

template <class T, eSignalDomain Domain>
auto BasicSignalView<T, Domain>::data() const -> pointer {
	return m_container.data();
}

template <class T, eSignalDomain Domain>
typename BasicSignalView<T, Domain>::size_type BasicSignalView<T, Domain>::size() const {
	return m_container.size();
}

template <class T, eSignalDomain Domain>
typename BasicSignalView<T, Domain>::size_type BasicSignalView<T, Domain>::size_bytes() const {
	return m_container.size_bytes();
}

template <class T, eSignalDomain Domain>
bool BasicSignalView<T, Domain>::empty() const {
	return m_container.empty();
}

template <class T, eSignalDomain Domain>
BasicSignalView<T, Domain> BasicSignalView<T, Domain>::first(size_type n) {
	return BasicSignalView(m_container.first(n));
}

template <class T, eSignalDomain Domain>
BasicSignalView<T, Domain> BasicSignalView<T, Domain>::last(size_type n) {
	return BasicSignalView(m_container.last(n));
}

template <class T, eSignalDomain Domain>
BasicSignalView<T, Domain> BasicSignalView<T, Domain>::subsignal(size_type offset) const {
	return BasicSignalView(m_container.subspan(offset));
}

template <class T, eSignalDomain Domain>
BasicSignalView<T, Domain> BasicSignalView<T, Domain>::subsignal(size_type offset, size_type count) const {
	return BasicSignalView(m_container.subspan(offset, count));
}

template <class T, eSignalDomain Domain>
auto BasicSignalView<T, Domain>::begin() const -> iterator {
	return m_container.begin();
}

template <class T, eSignalDomain Domain>
auto BasicSignalView<T, Domain>::end() const -> iterator {
	return m_container.end();
}

template <class T, eSignalDomain Domain>
auto BasicSignalView<T, Domain>::rbegin() const -> reverse_iterator {
	return m_container.rbegin();
}

template <class T, eSignalDomain Domain>
auto BasicSignalView<T, Domain>::rend() const -> reverse_iterator {
	return m_container.rend();
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


} // namespace dspbb

#include "Arithmetic.hpp"