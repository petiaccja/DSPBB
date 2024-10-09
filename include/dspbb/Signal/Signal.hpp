#pragma once

#include "Definitions.hpp"

#include <cassert>
#include <vector>


namespace dspbb {


template <class T, eSignalDomain Domain>
class BasicSignal {
	using container = std::vector<T>;

public:
	using value_type = typename container::value_type;
	using pointer = typename container::pointer;
	using const_pointer = typename container::const_pointer;
	using reference = typename container::reference;
	using const_reference = typename container::const_reference;
	using size_type = typename container::size_type;

	using iterator = typename container::iterator;
	using const_iterator = typename container::const_iterator;
	using reverse_iterator = typename container::reverse_iterator;
	using const_reverse_iterator = typename container::const_reverse_iterator;

public:
	BasicSignal() = default;
	BasicSignal(const BasicSignal&) = default;
	BasicSignal(BasicSignal&&) noexcept = default;
	BasicSignal& operator=(const BasicSignal&) = default;
	BasicSignal& operator=(BasicSignal&&) noexcept = default;

	explicit BasicSignal(container c) : m_container(std::move(c)) {}

	explicit BasicSignal(size_type count, const T& value = {});
	BasicSignal(std::initializer_list<T> ilist);
	template <class U>
	explicit BasicSignal(const BasicSignal<U, Domain>& other);
	BasicSignal(size_type count, const T* data);
	template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter>()), T>, int> = 0>
	BasicSignal(Iter first, Iter last) : m_container(first, last) {}

	template <class U>
	BasicSignal& operator=(const BasicSignal<U, Domain>&);

	reference operator[](size_t index);
	const_reference operator[](size_t index) const;
	pointer data();
	const_pointer data() const;

	size_type size() const;
	bool empty() const;
	size_type capacity() const;
	void reserve(size_type capacity);
	void resize(size_type count);
	void resize(size_type count, const T& value);

	void clear();
	void append(const BasicSignal& signal);
	void prepend(const BasicSignal& signal);
	void push_back(const T& value);
	BasicSignal extract_front(size_t count);
	BasicSignal extract_back(size_t count);
	void insert(size_type where, const BasicSignal& signal);
	void insert(const_iterator where, const BasicSignal& signal);
	template <class Iter>
	void insert(const_iterator where, Iter first, Iter last);
	void erase(const_iterator where);
	void erase(const_iterator first, const_iterator last);

	iterator begin();
	const_iterator begin() const;
	const_iterator cbegin() const;
	iterator end();
	const_iterator end() const;
	const_iterator cend() const;
	reverse_iterator rbegin();
	const_reverse_iterator rbegin() const;
	const_reverse_iterator crbegin() const;
	reverse_iterator rend();
	const_reverse_iterator rend() const;
	const_reverse_iterator crend() const;

private:
	container m_container;
};


//------------------------------------------------------------------------------
// Real signal
//------------------------------------------------------------------------------

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain>::BasicSignal(size_type count, const T& value) : m_container(count, value) {}

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain>::BasicSignal(std::initializer_list<T> ilist) : m_container(ilist) {}

template <class T, eSignalDomain Domain>
template <class U>
BasicSignal<T, Domain>::BasicSignal(const BasicSignal<U, Domain>& other) : m_container(other.begin(), other.end()) {
}

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain>::BasicSignal(size_type count, const T* data)
	: m_container(data, data + count) {}

template <class T, eSignalDomain Domain>
template <class U>
BasicSignal<T, Domain>& BasicSignal<T, Domain>::operator=(const BasicSignal<U, Domain>& other) {
	m_container.assign(other.begin(), other.end());
	return *this;
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::reference BasicSignal<T, Domain>::operator[](size_t index) {
	return m_container[index];
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::const_reference BasicSignal<T, Domain>::operator[](size_t index) const {
	return m_container[index];
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::pointer BasicSignal<T, Domain>::data() {
	return m_container.data();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::const_pointer BasicSignal<T, Domain>::data() const {
	return m_container.data();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::size_type BasicSignal<T, Domain>::size() const {
	return m_container.size();
}

template <class T, eSignalDomain Domain>
bool BasicSignal<T, Domain>::empty() const {
	return m_container.empty();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::size_type BasicSignal<T, Domain>::capacity() const {
	return m_container.capacity();
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::reserve(size_type capacity) {
	m_container.reserve(capacity);
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::resize(size_type count) {
	m_container.resize(count);
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::resize(size_type count, const T& value) {
	m_container.resize(count, value);
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::clear() {
	m_container.clear();
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::append(const BasicSignal& signal) {
	m_container.insert(m_container.end(), signal.begin(), signal.end());
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::prepend(const BasicSignal& signal) {
	m_container.insert(m_container.begin(), signal.begin(), signal.end());
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::push_back(const T& value) {
	m_container.push_back(value);
}

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain> BasicSignal<T, Domain>::extract_front(size_t count) {
	assert(count <= size());
	BasicSignal part{ count, data() };
	erase(begin(), begin() + count);
	return part;
}

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain> BasicSignal<T, Domain>::extract_back(size_t count) {
	assert(count <= size());
	BasicSignal part{ count, data() - count + size() };
	erase(end() - count, end());
	return part;
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::insert(size_type where, const BasicSignal& signal) {
	m_container.insert(m_container.begin() + where, signal.begin(), signal.end());
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::insert(const_iterator where, const BasicSignal& signal) {
	m_container.insert(where, signal.begin(), signal.end());
}

template <class T, eSignalDomain Domain>
template <class Iter>
void BasicSignal<T, Domain>::insert(const_iterator where, Iter first, Iter last) {
	m_container.insert(where, first, last);
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::erase(const_iterator where) {
	m_container.erase(where);
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::erase(const_iterator first, const_iterator last) {
	m_container.erase(first, last);
}

template <class T, eSignalDomain Domain>
auto BasicSignal<T, Domain>::begin() -> iterator {
	return m_container.begin();
}

template <class T, eSignalDomain Domain>
auto BasicSignal<T, Domain>::begin() const -> const_iterator {
	return m_container.begin();
}

template <class T, eSignalDomain Domain>
auto BasicSignal<T, Domain>::cbegin() const -> const_iterator {
	return m_container.cbegin();
}

template <class T, eSignalDomain Domain>
auto BasicSignal<T, Domain>::end() -> iterator {
	return m_container.end();
}

template <class T, eSignalDomain Domain>
auto BasicSignal<T, Domain>::end() const -> const_iterator {
	return m_container.end();
}

template <class T, eSignalDomain Domain>
auto BasicSignal<T, Domain>::cend() const -> const_iterator {
	return m_container.cend();
}

template <class T, eSignalDomain Domain>
auto BasicSignal<T, Domain>::rbegin() -> reverse_iterator {
	return m_container.rbegin();
}

template <class T, eSignalDomain Domain>
auto BasicSignal<T, Domain>::rbegin() const -> const_reverse_iterator {
	return m_container.rbegin();
}

template <class T, eSignalDomain Domain>
auto BasicSignal<T, Domain>::crbegin() const -> const_reverse_iterator {
	return m_container.crbegin();
}

template <class T, eSignalDomain Domain>
auto BasicSignal<T, Domain>::rend() -> reverse_iterator {
	return m_container.rend();
}

template <class T, eSignalDomain Domain>
auto BasicSignal<T, Domain>::rend() const -> const_reverse_iterator {
	return m_container.rend();
}

template <class T, eSignalDomain Domain>
auto BasicSignal<T, Domain>::crend() const -> const_reverse_iterator {
	return m_container.crend();
}

//------------------------------------------------------------------------------
// Helper types
//------------------------------------------------------------------------------

template <class T>
using Signal = BasicSignal<T, eSignalDomain::TIME>;
template <class T>
using Spectrum = BasicSignal<T, eSignalDomain::FREQUENCY>;
template <class T>
using Cepstrum = BasicSignal<T, eSignalDomain::QUEFRENCY>;


} // namespace dspbb


#include "Arithmetic.hpp"