#pragma once

#include "Signal.hpp"

#include "../Utility/TypeTraits.hpp"
#include "../Utility/TemplateUtil.hpp"


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
	SignalView(SignalT& signal, Q = 0) : SignalView(signal.begin(), signal.end()) {}
	template <class Q = std::enable_if_t<std::is_const<T>::value, int>>
	SignalView(const SignalT& signal, Q = 0) : SignalView(signal.begin(), signal.end()) {}
	SignalView(iterator first, iterator last) : first(first), last(last) {}
	SignalView(iterator first, size_t size) : first(first), last(first + size) {}

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

	decltype(auto) operator[](typename SignalT::size_type index) { return first[index]; }
	decltype(auto) operator[](typename SignalT::size_type index) const { return first[index]; }

	size_type Size() const { return size_type(last - first); }
	size_type Length() const { return Size(); }
	bool Empty() const { return first == last; }

	SignalView Subspan(size_type offset) const {
		assert(offset <= this->Size());
		return { this->first + offset, this->last };
	}
	SignalView Subspan(size_type offset, size_type count) const {
		assert(offset <= this->Size());
		assert(offset + count <= this->Size());
		return { this->first + offset, this->first + offset + count };
	}

	template <class U = std::enable_if_t<!std::is_const_v<T>, const T>>
	operator SignalView<U, Domain>() const {
		return SignalView<U, Domain>{ this->begin(), this->end() };
	}
	
	T* Data() { return std::addressof(*this->first); }
	const T* Data() const { return std::addressof(*this->first); }

protected:
	iterator first, last;
};

template <class T, eSignalDomain Domain>
SignalView<T, Domain> AsSpan(Signal<T, Domain>& signal) {
	return { signal.begin(), signal.end() };
}

template <class T, eSignalDomain Domain>
SignalView<const T, Domain> AsSpan(const Signal<T, Domain>& signal) {
	return { signal.begin(), signal.end() };
}

template <class T, eSignalDomain Domain>
SignalView<const T, Domain> AsConstSpan(const Signal<T, Domain>& signal) {
	return { signal.begin(), signal.end() };
}



} // namespace dspbb