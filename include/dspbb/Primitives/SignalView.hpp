#pragma once

#include "../Utility/TypeTraits.hpp"
#include "Signal.hpp"


namespace dspbb {


template <class T, eSignalDomain Domain>
class SignalView {
public:
	using SignalT = Signal<std::remove_const_t<T>, Domain>;

	using iterator = std::conditional_t<!std::is_const_v<T>, typename SignalT::iterator, typename SignalT::const_iterator>;
	using const_iterator = typename SignalT::const_iterator;
	using reverse_iterator = std::conditional_t<!std::is_const_v<T>, typename SignalT::reverse_iterator, typename SignalT::const_reverse_iterator>;
	using const_reverse_iterator = typename SignalT::const_reverse_iterator;
	using size_type = typename SignalT::size_type;

public:
	SignalView() = default;

	template <class Q = std::enable_if_t<!std::is_const<T>::value, int>>
	SignalView(SignalT& signal, Q = 0);

	template <class Q = std::enable_if_t<std::is_const<T>::value, int>>
	SignalView(const SignalT& signal, Q = 0);
	SignalView(iterator first, iterator last);
	SignalView(iterator first, size_t size);

	iterator begin() { return first; }
	const_iterator begin() const { return first; }
	const_iterator cbegin() const { return first; }
	iterator end() { return last; }
	const_iterator end() const { return last; }
	const_iterator cend() const { return last; }
	reverse_iterator rbegin() { return reverse_iterator{ first }; }
	const_reverse_iterator rbegin() const { return const_reverse_iterator{ first }; }
	const_reverse_iterator crbegin() const { return const_reverse_iterator{ first }; }
	reverse_iterator rend() { return reverse_iterator{ last }; }
	const_reverse_iterator rend() const { return const_reverse_iterator{ last }; }
	const_reverse_iterator crend() const { return const_reverse_iterator{ last }; }

	T& Front();
	const T& Front() const;
	T& Back();
	const T& Back() const;
	T& operator[](typename SignalT::size_type index);
	const T& operator[](typename SignalT::size_type index) const;
	T* Data();
	const T* Data() const;


	size_type Size() const;
	size_type Length() const;
	size_type SizeBytes() const;
	bool Empty() const;

	SignalView First(size_type n);
	SignalView Last(size_type n);
	SignalView SubSignal(size_type offset) const;
	SignalView SubSignal(size_type offset, size_type count) const;

	template <class U = std::enable_if_t<!std::is_const_v<T>, const T>>
	operator SignalView<U, Domain>() const;

	SignalView& operator+=(const SignalView& rhs);
	SignalView& operator-=(const SignalView& rhs);
	SignalView& operator*=(const SignalView& rhs);
	SignalView& operator/=(const SignalView& rhs);

	SignalView& operator+=(const T& scalar);
	SignalView& operator-=(const T& scalar);
	SignalView& operator*=(const T& scalar);
	SignalView& operator/=(const T& scalar);

	friend SignalT operator+=(SignalT lhs, const SignalView& rhs) {
		assert(lhs.Size() == rhs.Size());
		for (size_t i = 0; i < lhs.Size(); ++i) {
			lhs[i] += rhs[i];
		}
		return lhs;
	}
	friend SignalT operator-=(SignalT lhs, const SignalView& rhs) {
		assert(lhs.Size() == rhs.Size());
		for (size_t i = 0; i < lhs.Size(); ++i) {
			lhs[i] -= rhs[i];
		}
		return lhs;
	}
	friend SignalT operator*=(SignalT lhs, const SignalView& rhs) {
		assert(lhs.Size() == rhs.Size());
		for (size_t i = 0; i < lhs.Size(); ++i) {
			lhs[i] *= rhs[i];
		}
		return lhs;
	}
	friend SignalT operator/=(SignalT lhs, const SignalView& rhs) {
		assert(lhs.Size() == rhs.Size());
		for (size_t i = 0; i < lhs.Size(); ++i) {
			lhs[i] /= rhs[i];
		}
		return lhs;
	}

	SignalView& operator+=(const SignalT& rhs) { return *this += SignalView(rhs); }
	SignalView& operator-=(const SignalT& rhs) { return *this -= SignalView(rhs); }
	SignalView& operator*=(const SignalT& rhs) { return *this *= SignalView(rhs); }
	SignalView& operator/=(const SignalT& rhs) { return *this /= SignalView(rhs); }

	SignalT operator+(const SignalView& rhs) const { return SignalT(begin(), end()) += rhs; }
	SignalT operator-(const SignalView& rhs) const { return SignalT(begin(), end()) -= rhs; }
	SignalT operator*(const SignalView& rhs) const { return SignalT(begin(), end()) *= rhs; }
	SignalT operator/(const SignalView& rhs) const { return SignalT(begin(), end()) /= rhs; }

	SignalT operator+(const T& scalar) const { return SignalT(begin(), end()) += scalar; }
	SignalT operator-(const T& scalar) const { return SignalT(begin(), end()) -= scalar; }
	SignalT operator*(const T& scalar) const { return SignalT(begin(), end()) *= scalar; }
	SignalT operator/(const T& scalar) const { return SignalT(begin(), end()) /= scalar; }

	friend SignalT operator+(const T& scalar, const SignalView& signal) { return scalar + SignalT(signal.begin(), signal.end()); }
	friend SignalT operator-(const T& scalar, const SignalView& signal) { return scalar - SignalT(signal.begin(), signal.end()); }
	friend SignalT operator*(const T& scalar, const SignalView& signal) { return scalar * SignalT(signal.begin(), signal.end()); }
	friend SignalT operator/(const T& scalar, const SignalView& signal) { return scalar / SignalT(signal.begin(), signal.end()); }

	friend SignalT operator+(SignalT lhs, const SignalView& rhs) { return lhs + SignalT(rhs.begin(), rhs.end()); }
	friend SignalT operator-(SignalT lhs, const SignalView& rhs) { return lhs - SignalT(rhs.begin(), rhs.end()); }
	friend SignalT operator*(SignalT lhs, const SignalView& rhs) { return lhs * SignalT(rhs.begin(), rhs.end()); }
	friend SignalT operator/(SignalT lhs, const SignalView& rhs) { return lhs / SignalT(rhs.begin(), rhs.end()); }



protected:
	iterator first, last;
};


template <class T, eSignalDomain Domain>
template <class Q>
SignalView<T, Domain>::SignalView(SignalT& signal, Q)
	: SignalView(signal.begin(), signal.end()) {
}

template <class T, eSignalDomain Domain>
template <class Q>
SignalView<T, Domain>::SignalView(const SignalT& signal, Q)
	: SignalView(signal.begin(), signal.end()) {
}

template <class T, eSignalDomain Domain>
SignalView<T, Domain>::SignalView(iterator first, iterator last)
	: first(first),
	  last(last) {
}

template <class T, eSignalDomain Domain>
SignalView<T, Domain>::SignalView(iterator first, size_t size)
	: first(first),
	  last(first + size) {
}

template <class T, eSignalDomain Domain>
T& SignalView<T, Domain>::Front() { return *first; }

template <class T, eSignalDomain Domain>
const T& SignalView<T, Domain>::Front() const { return *first; }

template <class T, eSignalDomain Domain>
T& SignalView<T, Domain>::Back() { return *(last - 1); }

template <class T, eSignalDomain Domain>
const T& SignalView<T, Domain>::Back() const { return *(last - 1); }

template <class T, eSignalDomain Domain>
T& SignalView<T, Domain>::operator[](typename SignalT::size_type index) { return first[index]; }

template <class T, eSignalDomain Domain>
const T& SignalView<T, Domain>::operator[](typename SignalT::size_type index) const { return first[index]; }

template <class T, eSignalDomain Domain>
T* SignalView<T, Domain>::Data() { return std::addressof(*first); }

template <class T, eSignalDomain Domain>
const T* SignalView<T, Domain>::Data() const { return std::addressof(*first); }

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

template <class T, eSignalDomain Domain>
template <class U>
SignalView<T, Domain>::operator SignalView<U, Domain>() const {
	return SignalView<U, Domain>{ this->begin(), this->end() };
}

template <class T, eSignalDomain Domain>
SignalView<T, Domain>& SignalView<T, Domain>::operator+=(const SignalView& rhs) {
	assert((*this).Size() == rhs.Size());
	for (size_t i = 0; i < (*this).Size(); ++i) {
		(*this)[i] += rhs[i];
	}
	return *this;
}

template <class T, eSignalDomain Domain>
SignalView<T, Domain>& SignalView<T, Domain>::operator-=(const SignalView& rhs) {
	assert((*this).Size() == rhs.Size());
	for (size_t i = 0; i < (*this).Size(); ++i) {
		(*this)[i] -= rhs[i];
	}
	return *this;
}

template <class T, eSignalDomain Domain>
SignalView<T, Domain>& SignalView<T, Domain>::operator*=(const SignalView& rhs) {
	assert((*this).Size() == rhs.Size());
	for (size_t i = 0; i < (*this).Size(); ++i) {
		(*this)[i] *= rhs[i];
	}
	return *this;
}

template <class T, eSignalDomain Domain>
SignalView<T, Domain>& SignalView<T, Domain>::operator/=(const SignalView& rhs) {
	assert((*this).Size() == rhs.Size());
	for (size_t i = 0; i < (*this).Size(); ++i) {
		(*this)[i] /= rhs[i];
	}
	return *this;
}

template <class T, eSignalDomain Domain>
SignalView<T, Domain>& SignalView<T, Domain>::operator+=(const T& scalar) {
	for (auto& v : *this) {
		v += scalar;
	}
	return *this;
}

template <class T, eSignalDomain Domain>
SignalView<T, Domain>& SignalView<T, Domain>::operator-=(const T& scalar) {
	for (auto& v : *this) {
		v -= scalar;
	}
	return *this;
}

template <class T, eSignalDomain Domain>
SignalView<T, Domain>& SignalView<T, Domain>::operator*=(const T& scalar) {
	for (auto& v : *this) {
		v *= scalar;
	}
	return *this;
}

template <class T, eSignalDomain Domain>
SignalView<T, Domain>& SignalView<T, Domain>::operator/=(const T& scalar) {
	for (auto& v : *this) {
		v /= scalar;
	}
	return *this;
}

// Helpers
template <class T, eSignalDomain Domain>
SignalView<T, Domain> AsView(Signal<T, Domain>& signal) {
	return { signal.begin(), signal.end() };
}

template <class T, eSignalDomain Domain>
SignalView<const T, Domain> AsView(const Signal<T, Domain>& signal) {
	return { signal.begin(), signal.end() };
}

template <class T, eSignalDomain Domain>
SignalView<const T, Domain> AsConstView(const Signal<T, Domain>& signal) {
	return { signal.begin(), signal.end() };
}

template <class T, eSignalDomain Domain>
SignalView<T, Domain> AsView(SignalView<T, Domain> view) {
	return view;
}

template <class T, eSignalDomain Domain>
SignalView<const T, Domain> AsView(SignalView<const T, Domain> view) {
	return view;
}

template <class T, eSignalDomain Domain>
SignalView<const T, Domain> AsConstView(SignalView<T, Domain> view) {
	return view;
}

template <class T, eSignalDomain Domain>
SignalView<const T, Domain> AsConstView(SignalView<const T, Domain> view) {
	return view;
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