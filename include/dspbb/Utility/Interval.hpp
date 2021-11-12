#pragma once

#include <algorithm>
#include <optional>

namespace dspbb {

template <class T>
struct Interval {
	Interval() = default;
	Interval(T first, T last) : first(first), last(last) {}
	T first;
	T last;
	T Size() const { return last - first; }
};

template <class T>
Interval<T> operator+(const Interval<T>& lhs, T rhs) {
	return { lhs.first + rhs, lhs.last + rhs };
}

template <class T>
Interval<T>& operator+=(Interval<T>& lhs, T rhs) {
	return lhs = lhs + rhs;
}

template <class T>
Interval<T> operator-(const Interval<T>& lhs, T rhs) {
	return { lhs.first - rhs, lhs.last - rhs };
}

template <class T>
Interval<T>& operator-=(Interval<T>& lhs, T rhs) {
	return lhs = lhs - rhs;
}


template <class T>
bool operator==(const Interval<T>& lhs, const Interval<T>& rhs) {
	return lhs.first == rhs.first && lhs.last && rhs.last;
}

template <class T>
bool operator!=(const Interval<T>& lhs, const Interval<T>& rhs) {
	return !(lhs == rhs);
}

template <class T>
bool IsDisjoint(const Interval<T>& lhs, const Interval<T>& rhs) {
	const T first = std::max(lhs.first, rhs.first);
	const T last = std::min(lhs.last, rhs.last);
	return !(first < last);
}

template <class T>
Interval<T> Intersection(const Interval<T>& lhs, const Interval<T>& rhs) {
	const T first = std::max(lhs.first, rhs.first);
	const T last = std::min(lhs.last, rhs.last);
	if (first < last) {
		return { first, last };
	}
	return { T(0), T(0) };
}

template <class T>
Interval<T> EncompassingUnion(const Interval<T>& lhs, const Interval<T>& rhs) {
	const T first = std::min(lhs.first, rhs.first);
	const T last = std::max(lhs.last, rhs.last);
	return { first, last };
}


template <class T>
std::optional<Interval<T>> Union(const Interval<T>& lhs, const Interval<T>& rhs) {
	if (!IsDisjoint(lhs, rhs)) {
		return { EncompassingUnion(lhs, rhs) };
	}
	return {};
}


} // namespace dspbb