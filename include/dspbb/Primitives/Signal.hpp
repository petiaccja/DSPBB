#pragma once

#include <cassert>
#include <complex>
#include <vector>


namespace dspbb {


enum class eSignalDomain {
	TIME,
	FREQUENCY,
	QUEFRENCY,
	DOMAINLESS,
};
static constexpr auto TIME_DOMAIN = eSignalDomain::TIME;
static constexpr auto FREQUENCY_DOMAIN = eSignalDomain::FREQUENCY;
static constexpr auto QUEFRENCY_DOMAIN = eSignalDomain::QUEFRENCY;
static constexpr auto DOMAINLESS = eSignalDomain::DOMAINLESS;


template <class T, eSignalDomain Domain>
class BasicSignal {
	template <class U, eSignalDomain DomainB>
	friend class BasicSignal;
	using storage_type = std::vector<T>;

public:
	using value_type = T;
	using pointer = T*;
	using const_pointer = const T*;
	using reference = value_type&;
	using const_reference = const value_type&;
	using size_type = std::size_t;

	using iterator = typename storage_type::iterator;
	using const_iterator = typename storage_type::const_iterator;
	using reverse_iterator = typename storage_type::reverse_iterator;
	using const_reverse_iterator = typename storage_type::const_reverse_iterator;

public:
	BasicSignal() = default;
	explicit BasicSignal(size_type count, const T& value = {});
	BasicSignal(const BasicSignal&) = default;
	BasicSignal(BasicSignal&&) noexcept = default;
	BasicSignal(std::initializer_list<T> ilist);
	template <class U>
	explicit BasicSignal(const BasicSignal<U, Domain>& other);
	BasicSignal(size_type count, const T* data);
	template <class Iter, std::enable_if_t<std::is_convertible_v<decltype(*std::declval<Iter>()), T>, int> = 0>
	BasicSignal(Iter first, Iter last) : m_samples(first, last) {}

	BasicSignal& operator=(const BasicSignal&) = default;
	BasicSignal& operator=(BasicSignal&&) noexcept = default;
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
	storage_type m_samples;
};


//------------------------------------------------------------------------------
// Real signal
//------------------------------------------------------------------------------

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain>::BasicSignal(size_type count, const T& value) : m_samples(count, value) {}

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain>::BasicSignal(std::initializer_list<T> ilist) : m_samples(ilist) {}

template <class T, eSignalDomain Domain>
template <class U>
BasicSignal<T, Domain>::BasicSignal(const BasicSignal<U, Domain>& other) : m_samples(other.begin(), other.end()) {
}

template <class T, eSignalDomain Domain>
BasicSignal<T, Domain>::BasicSignal(size_type count, const T* data)
	: m_samples(data, data + count) {}

template <class T, eSignalDomain Domain>
template <class U>
BasicSignal<T, Domain>& BasicSignal<T, Domain>::operator=(const BasicSignal<U, Domain>& other) {
	m_samples.assign(other.begin(), other.end());
	return *this;
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::reference BasicSignal<T, Domain>::operator[](size_t index) {
	return m_samples[index];
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::const_reference BasicSignal<T, Domain>::operator[](size_t index) const {
	return m_samples[index];
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::pointer BasicSignal<T, Domain>::data() {
	return m_samples.data();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::const_pointer BasicSignal<T, Domain>::data() const {
	return m_samples.data();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::size_type BasicSignal<T, Domain>::size() const {
	return m_samples.size();
}

template <class T, eSignalDomain Domain>
bool BasicSignal<T, Domain>::empty() const {
	return m_samples.empty();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::size_type BasicSignal<T, Domain>::capacity() const {
	return m_samples.capacity();
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::reserve(size_type capacity) {
	m_samples.reserve(capacity);
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::resize(size_type count) {
	m_samples.resize(count);
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::resize(size_type count, const T& value) {
	m_samples.resize(count, value);
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::clear() {
	m_samples.clear();
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::append(const BasicSignal& signal) {
	m_samples.insert(m_samples.end(), signal.begin(), signal.end());
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::prepend(const BasicSignal& signal) {
	m_samples.insert(m_samples.begin(), signal.begin(), signal.end());
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::push_back(const T& value) {
	m_samples.push_back(value);
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
	m_samples.insert(m_samples.begin() + where, signal.begin(), signal.end());
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::insert(const_iterator where, const BasicSignal& signal) {
	m_samples.insert(where, signal.begin(), signal.end());
}

template <class T, eSignalDomain Domain>
template <class Iter>
void BasicSignal<T, Domain>::insert(const_iterator where, Iter first, Iter last) {
	m_samples.insert(where, first, last);
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::erase(const_iterator where) {
	m_samples.erase(where);
}

template <class T, eSignalDomain Domain>
void BasicSignal<T, Domain>::erase(const_iterator first, const_iterator last) {
	m_samples.erase(first, last);
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::iterator BasicSignal<T, Domain>::begin() {
	return m_samples.begin();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::const_iterator BasicSignal<T, Domain>::begin() const {
	return m_samples.begin();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::const_iterator BasicSignal<T, Domain>::cbegin() const {
	return m_samples.cbegin();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::iterator BasicSignal<T, Domain>::end() {
	return m_samples.end();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::const_iterator BasicSignal<T, Domain>::end() const {
	return m_samples.end();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::const_iterator BasicSignal<T, Domain>::cend() const {
	return m_samples.cend();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::reverse_iterator BasicSignal<T, Domain>::rbegin() {
	return m_samples.rbegin();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::const_reverse_iterator BasicSignal<T, Domain>::rbegin() const {
	return m_samples.rbegin();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::const_reverse_iterator BasicSignal<T, Domain>::crbegin() const {
	return m_samples.crbegin();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::reverse_iterator BasicSignal<T, Domain>::rend() {
	return m_samples.rend();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::const_reverse_iterator BasicSignal<T, Domain>::rend() const {
	return m_samples.rend();
}

template <class T, eSignalDomain Domain>
typename BasicSignal<T, Domain>::const_reverse_iterator BasicSignal<T, Domain>::crend() const {
	return m_samples.crend();
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

using SignalF = Signal<float>;
using SignalCF = Signal<std::complex<float>>;

using SpectrumCF = BasicSignal<std::complex<float>, eSignalDomain::FREQUENCY>;
using SpectrumF = BasicSignal<float, eSignalDomain::FREQUENCY>;

using CepstrumCF = BasicSignal<std::complex<float>, eSignalDomain::QUEFRENCY>;
using CepstrumF = BasicSignal<float, eSignalDomain::QUEFRENCY>;


} // namespace dspbb


#include "SignalArithmetic.hpp"