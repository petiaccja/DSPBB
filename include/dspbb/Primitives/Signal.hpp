#pragma once

#include <algorithm>
#include <complex>
#include <vector>
#include <cassert>
#include <stdexcept>


namespace enl {


enum class eSignalDomain {
	TIME,
	FREQUENCY,
	QUEFRENCY,
	DOMAINLESS,
};
static constexpr auto TIME_DOMAIN = eSignalDomain::TIME;
static constexpr auto FREQUENCY_DOMAIN = eSignalDomain::FREQUENCY;
static constexpr auto QUEFRENCY_DOMAIN = eSignalDomain::QUEFRENCY;


template <class T, eSignalDomain Domain>
class Signal {
	template <class U, eSignalDomain DomainB>
	friend class Signal;

public:
	using value_type = T;
	using pointer = T*;
	using const_pointer = const T*;
	using reference = value_type&;
	using const_reference = const value_type&;
	using size_type = std::size_t;

	using iterator = typename std::vector<T>::iterator;
	using const_iterator = typename std::vector<T>::const_iterator;
	using reverse_iterator = typename std::vector<T>::reverse_iterator;
	using const_reverse_iterator = typename std::vector<T>::const_reverse_iterator;

public:
	Signal() = default;
	explicit Signal(size_type count, const T& value = {});
	Signal(const Signal&) = default;
	Signal(Signal&&) = default;
	Signal(std::initializer_list<T> ilist);
	template <class U>
	Signal(const Signal<U, Domain>& other);
	Signal(size_type count, const T* data);
	template <class Iter>
	Signal(Iter first, Iter last) : m_samples(first, last) {}

	Signal& operator=(const Signal&) = default;
	Signal& operator=(Signal&&) = default;
	template <class U>
	Signal& operator=(const Signal<U, Domain>&);

	reference operator[](size_t index);
	const_reference operator[](size_t index) const;
	pointer Data();
	const_pointer Data() const;

	size_type Size() const;
	size_type Length() const;
	bool Empty() const;
	size_type Capacity() const;
	void Reserve(size_type capacity);
	void Resize(size_type count);
	void Resize(size_type count, const T& value);

	void Clear();
	void Append(const Signal& signal);
	void Prepend(const Signal& signal);
	void PushBack(const T& value);
	Signal ExtractFront(size_t count);
	Signal ExtractBack(size_t count);
	void Insert(size_type where, const Signal& signal);
	void Insert(const_iterator where, const Signal& signal);
	template <class Iter>
	void Insert(const_iterator where, Iter first, Iter last);
	void Erase(const_iterator where);
	void Erase(const_iterator first, const_iterator last);

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

	Signal& operator+=(const Signal& rhs);
	Signal& operator-=(const Signal& rhs);
	Signal& operator*=(const Signal& rhs);
	Signal& operator/=(const Signal& rhs);
	Signal& operator+=(const T& scalar);
	Signal& operator-=(const T& scalar);
	Signal& operator*=(const T& scalar);
	Signal& operator/=(const T& scalar);

	Signal operator+(const Signal& rhs) const { return Signal(*this) += rhs; }
	Signal operator-(const Signal& rhs) const { return Signal(*this) -= rhs; }
	Signal operator*(const Signal& rhs) const { return Signal(*this) *= rhs; }
	Signal operator/(const Signal& rhs) const { return Signal(*this) /= rhs; }
	Signal operator+(const T& scalar) const { return Signal(*this) += scalar; }
	Signal operator-(const T& scalar) const { return Signal(*this) -= scalar; }
	Signal operator*(const T& scalar) const { return Signal(*this) *= scalar; }
	Signal operator/(const T& scalar) const { return Signal(*this) /= scalar; }

	friend Signal operator+(const T& scalar, const Signal& signal) { return signal * scalar; }
	friend Signal operator-(const T& scalar, const Signal& signal) { return signal * scalar; }
	friend Signal operator*(const T& scalar, const Signal& signal) { return signal * scalar; }
	friend Signal operator/(const T& scalar, const Signal& signal) { return signal * scalar; }

private:
	std::vector<T> m_samples;
};



//------------------------------------------------------------------------------
// Real signal
//------------------------------------------------------------------------------

template <class T, eSignalDomain Domain>
Signal<T, Domain>::Signal(size_type count, const T& value) : m_samples(count, value) {}

template <class T, eSignalDomain Domain>
Signal<T, Domain>::Signal(std::initializer_list<T> ilist) : m_samples(ilist) {}

template <class T, eSignalDomain Domain>
template <class U>
Signal<T, Domain>::Signal(const Signal<U, Domain>& other) : m_samples(other.begin(), other.end()) {
}

template <class T, eSignalDomain Domain>
Signal<T, Domain>::Signal(size_type count, const T* data)
	: m_samples(data, data + count) {}

template <class T, eSignalDomain Domain>
template <class U>
Signal<T, Domain>& Signal<T, Domain>::operator=(const Signal<U, Domain>& other) {
	m_samples.assign(other.begin(), other.end());
	return *this;
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::reference Signal<T, Domain>::operator[](size_t index) {
	return m_samples[index];
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::const_reference Signal<T, Domain>::operator[](size_t index) const {
	return m_samples[index];
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::pointer Signal<T, Domain>::Data() {
	return m_samples.data();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::const_pointer Signal<T, Domain>::Data() const {
	return m_samples.data();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::size_type Signal<T, Domain>::Size() const {
	return m_samples.size();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::size_type Signal<T, Domain>::Length() const {
	return Size();
}

template <class T, eSignalDomain Domain>
bool Signal<T, Domain>::Empty() const {
	return m_samples.empty();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::size_type Signal<T, Domain>::Capacity() const {
	return m_samples.capacity();
}

template <class T, eSignalDomain Domain>
void Signal<T, Domain>::Reserve(size_type capacity) {
	m_samples.reserve(capacity);
}

template <class T, eSignalDomain Domain>
void Signal<T, Domain>::Resize(size_type count) {
	m_samples.resize(count);
}

template <class T, eSignalDomain Domain>
void Signal<T, Domain>::Resize(size_type count, const T& value) {
	m_samples.resize(count, value);
}

template <class T, eSignalDomain Domain>
void Signal<T, Domain>::Clear() {
	m_samples.clear();
}

template <class T, eSignalDomain Domain>
void Signal<T, Domain>::Append(const Signal& signal) {
	m_samples.insert(m_samples.end(), signal.begin(), signal.end());
}

template <class T, eSignalDomain Domain>
void Signal<T, Domain>::Prepend(const Signal& signal) {
	m_samples.insert(m_samples.begin(), signal.begin(), signal.end());
}

template <class T, eSignalDomain Domain>
void Signal<T, Domain>::PushBack(const T& value) {
	m_samples.push_back(value);
}

template <class T, eSignalDomain Domain>
Signal<T, Domain> Signal<T, Domain>::ExtractFront(size_t count) {
	assert(count <= Length());
	Signal part{ count, Data() };
	Erase(begin(), begin() + count);
	return part;
}

template <class T, eSignalDomain Domain>
Signal<T, Domain> Signal<T, Domain>::ExtractBack(size_t count) {
	assert(count <= Length());
	Signal part{ count, Data() - count + Length() };
	Erase(end() - count, end());
	return part;
}

template <class T, eSignalDomain Domain>
void Signal<T, Domain>::Insert(size_type where, const Signal& signal) {
	m_samples.insert(m_samples.begin() + where, signal.begin(), signal.end());
}

template <class T, eSignalDomain Domain>
void Signal<T, Domain>::Insert(const_iterator where, const Signal& signal) {
	m_samples.insert(where, signal.begin(), signal.end());
}

template <class T, eSignalDomain Domain>
template <class Iter>
void Signal<T, Domain>::Insert(const_iterator where, Iter first, Iter last) {
	m_samples.insert(where, first, last);
}

template <class T, eSignalDomain Domain>
void Signal<T, Domain>::Erase(const_iterator where) {
	m_samples.erase(where);
}

template <class T, eSignalDomain Domain>
void Signal<T, Domain>::Erase(const_iterator first, const_iterator last) {
	m_samples.erase(first, last);
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::iterator Signal<T, Domain>::begin() {
	return m_samples.begin();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::const_iterator Signal<T, Domain>::begin() const {
	return m_samples.begin();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::const_iterator Signal<T, Domain>::cbegin() const {
	return m_samples.cbegin();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::iterator Signal<T, Domain>::end() {
	return m_samples.end();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::const_iterator Signal<T, Domain>::end() const {
	return m_samples.end();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::const_iterator Signal<T, Domain>::cend() const {
	return m_samples.cend();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::reverse_iterator Signal<T, Domain>::rbegin() {
	return m_samples.rbegin();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::const_reverse_iterator Signal<T, Domain>::rbegin() const {
	return m_samples.rbegin();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::const_reverse_iterator Signal<T, Domain>::crbegin() const {
	return m_samples.crbegin();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::reverse_iterator Signal<T, Domain>::rend() {
	return m_samples.rend();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::const_reverse_iterator Signal<T, Domain>::rend() const {
	return m_samples.rend();
}

template <class T, eSignalDomain Domain>
typename Signal<T, Domain>::const_reverse_iterator Signal<T, Domain>::crend() const {
	return m_samples.crend();
}

template <class T, eSignalDomain Domain>
Signal<T, Domain>& Signal<T, Domain>::operator+=(const Signal& rhs) {
	if (Length() != rhs.Length()) {
		throw std::logic_error("Signals must be exactly the same lenth.");
	}
	for (size_t i = 0; i < Length(); ++i) {
		(*this)[i] += rhs[i];
	}
	return *this;
}

template <class T, eSignalDomain Domain>
Signal<T, Domain>& Signal<T, Domain>::operator-=(const Signal& rhs) {
	if (Length() != rhs.Length()) {
		throw std::logic_error("Signals must be exactly the same lenth.");
	}
	for (size_t i =0; i<Length(); ++i) {
		(*this)[i] -= rhs[i];
	}
	return *this;
}

template <class T, eSignalDomain Domain>
Signal<T, Domain>& Signal<T, Domain>::operator*=(const Signal& rhs) {
	if (Length() != rhs.Length()) {
		throw std::logic_error("Signals must be exactly the same lenth.");
	}
	for (size_t i = 0; i < Length(); ++i) {
		(*this)[i] *= rhs[i];
	}
	return *this;
}

template <class T, eSignalDomain Domain>
Signal<T, Domain>& Signal<T, Domain>::operator/=(const Signal& rhs) {
	if (Length() != rhs.Length()) {
		throw std::logic_error("Signals must be exactly the same lenth.");
	}
	for (size_t i = 0; i < Length(); ++i) {
		(*this)[i] /= rhs[i];
	}
	return *this;
}

template <class T, eSignalDomain Domain>
Signal<T, Domain>& Signal<T, Domain>::operator+=(const T& scalar) {
	for (size_t i = 0; i < Length(); ++i) {
		(*this)[i] += scalar;
	}
	return *this;
}

template <class T, eSignalDomain Domain>
Signal<T, Domain>& Signal<T, Domain>::operator-=(const T& scalar) {
	for (size_t i = 0; i < Length(); ++i) {
		(*this)[i] -= scalar;
	}
	return *this;
}

template <class T, eSignalDomain Domain>
Signal<T, Domain>& Signal<T, Domain>::operator*=(const T& scalar) {
	for (size_t i = 0; i < Length(); ++i) {
		(*this)[i] *= scalar;
	}
	return *this;
}

template <class T, eSignalDomain Domain>
Signal<T, Domain>& Signal<T, Domain>::operator/=(const T& scalar) {
	for (size_t i = 0; i < Length(); ++i) {
		(*this)[i] /= scalar;
	}
	return *this;
}

//------------------------------------------------------------------------------
// Helper types
//------------------------------------------------------------------------------

template <class T>
using TimeSignal = Signal<T, eSignalDomain::TIME>;
template <class T>
using Spectrum = Signal<T, eSignalDomain::FREQUENCY>;
template <class T>
using Cepstrum = Signal<T, eSignalDomain::QUEFRENCY>;

using TimeSignalF = TimeSignal<float>;
using TimeSignalCF = TimeSignal<std::complex<float>>;

using SpectrumCF = Signal<std::complex<float>, eSignalDomain::FREQUENCY>;
using SpectrumF = Signal<float, eSignalDomain::FREQUENCY>;


} // namespace enl