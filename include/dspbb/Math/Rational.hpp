#pragma once

#include <numeric>
#include <tuple>
#include <type_traits>

namespace dspbb {

template <class T>
class Rational {
	static_assert(std::is_integral_v<T>);

public:
	using int_t = T;

public:
	constexpr Rational() noexcept = default;
	constexpr Rational(const Rational&) noexcept = default;
	constexpr explicit Rational(int_t value) noexcept;
	constexpr Rational(int_t numerator, int_t denominator) noexcept;
	Rational& operator=(const Rational&) = default;

	constexpr int_t Numerator() const;
	constexpr int_t Denominator() const;

	Rational& operator++() noexcept;
	Rational& operator--() noexcept;
	Rational operator++(int) noexcept;
	Rational operator--(int) noexcept;

	Rational& operator+=(int_t rhs) noexcept;
	Rational& operator-=(int_t rhs) noexcept;
	Rational& operator*=(int_t rhs) noexcept;
	Rational& operator/=(int_t rhs) noexcept;

	Rational& operator+=(const Rational& rhs) noexcept;
	Rational& operator-=(const Rational& rhs) noexcept;
	Rational& operator*=(const Rational& rhs) noexcept;
	Rational& operator/=(const Rational& rhs) noexcept;

	template <class FloatT, std::enable_if_t<std::is_floating_point_v<FloatT>, int> = 0>
	operator FloatT() const;

private:
	int_t num;
	int_t den;
};


template <class T>
constexpr Rational<T> operator+(const Rational<T>& lhs, T rhs) noexcept {
	return { lhs.Numerator() + lhs.Denominator() * rhs,
			 lhs.Denominator() };
}
template <class T>
constexpr Rational<T> operator-(const Rational<T>& lhs, T rhs) noexcept {
	return { lhs.Numerator() - lhs.Denominator() * rhs,
			 lhs.Denominator() };
}
template <class T>
constexpr Rational<T> operator*(const Rational<T>& lhs, T rhs) noexcept {
	const auto simplification = std::gcd(rhs, lhs.Denominator());
	return { lhs.Numerator() * (rhs / simplification),
			 lhs.Denominator() / simplification };
}
template <class T>
constexpr Rational<T> operator/(const Rational<T>& lhs, T rhs) noexcept {
	const auto simpl = std::gcd(rhs, lhs.Numerator());
	return { lhs.Numerator() / simpl,
			 lhs.Denominator() * (rhs / simpl) };
}


template <class T>
constexpr Rational<T> operator+(const Rational<T>& lhs, const Rational<T>& rhs) noexcept {
	const auto commonDen = std::lcm(lhs.Denominator(), rhs.Denominator());
	const auto num = lhs.Numerator() * (commonDen / lhs.Denominator()) + rhs.Numerator() * (commonDen / rhs.Denominator());
	const auto den = commonDen;
	const auto simpl = std::gcd(num, den);
	return { num / simpl, den / simpl };
}
template <class T>
constexpr Rational<T> operator-(const Rational<T>& lhs, const Rational<T>& rhs) noexcept {
	const auto commonDen = std::lcm(lhs.Denominator(), rhs.Denominator());
	const auto num = lhs.Numerator() * (commonDen / lhs.Denominator()) - rhs.Numerator() * (commonDen / rhs.Denominator());
	const auto den = commonDen;
	const auto simpl = std::gcd(num, den);
	return { num / simpl, den / simpl };
}
template <class T>
constexpr Rational<T> operator*(const Rational<T>& lhs, const Rational<T>& rhs) noexcept {
	const auto simpl1 = std::gcd(lhs.Numerator(), rhs.Denominator());
	const auto simpl2 = std::gcd(rhs.Numerator(), lhs.Denominator());
	return { (lhs.Numerator() / simpl1) * (rhs.Numerator() / simpl2),
			 (rhs.Denominator() / simpl1) * (lhs.Denominator() / simpl2) };
}
template <class T>
constexpr Rational<T> operator/(const Rational<T>& lhs, const Rational<T>& rhs) noexcept {
	return lhs * Rational<T>{ rhs.Denominator(), rhs.Numerator() };
}


template <class T>
constexpr Rational<T>::Rational(int_t value) noexcept
	: num(value), den(1) {}

template <class T>
constexpr Rational<T>::Rational(int_t numerator, int_t denominator) noexcept
	: num(numerator / std::gcd(numerator, denominator)), den(denominator / std::gcd(numerator, denominator)) {}

template <class T>
constexpr typename Rational<T>::int_t Rational<T>::Numerator() const {
	return num;
}

template <class T>
constexpr typename Rational<T>::int_t Rational<T>::Denominator() const {
	return den;
}

template <class T>
Rational<T>& Rational<T>::operator++() noexcept {
	num += den;
	return *this;
}

template <class T>
Rational<T>& Rational<T>::operator--() noexcept {
	num -= den;
	return *this;
}

template <class T>
Rational<T> Rational<T>::operator++(int) noexcept {
	const auto copy = *this;
	++*this;
	return copy;
}

template <class T>
Rational<T> Rational<T>::operator--(int) noexcept {
	const auto copy = *this;
	--*this;
	return copy;
}


template <class T>
Rational<T>& Rational<T>::operator+=(int_t rhs) noexcept {
	return *this = *this + rhs;
}

template <class T>
Rational<T>& Rational<T>::operator-=(int_t rhs) noexcept {
	return *this = *this - rhs;
}

template <class T>
Rational<T>& Rational<T>::operator*=(int_t rhs) noexcept {
	return *this = *this * rhs;
}

template <class T>
Rational<T>& Rational<T>::operator/=(int_t rhs) noexcept {
	return *this = *this / rhs;
}



template <class T>
Rational<T>& Rational<T>::operator+=(const Rational& rhs) noexcept {
	return *this = *this + rhs;
}

template <class T>
Rational<T>& Rational<T>::operator-=(const Rational& rhs) noexcept {
	return *this = *this - rhs;
}

template <class T>
Rational<T>& Rational<T>::operator*=(const Rational& rhs) noexcept {
	return *this = *this * rhs;
}

template <class T>
Rational<T>& Rational<T>::operator/=(const Rational& rhs) noexcept {
	return *this = *this / rhs;
}

template <class T>
template <class FloatT, std::enable_if_t<std::is_floating_point_v<FloatT>, int>>
Rational<T>::operator FloatT() const {
	return static_cast<FloatT>(num) / static_cast<FloatT>(den);
}


template <class T>
constexpr Rational<T> operator+(T lhs, const Rational<T>& rhs) noexcept {
	return rhs + lhs;
}
template <class T>
constexpr Rational<T> operator-(T lhs, const Rational<T>& rhs) noexcept {
	return Rational<T>{ lhs } - rhs;
}
template <class T>
constexpr Rational<T> operator*(T lhs, const Rational<T>& rhs) noexcept {
	return Rational<T>{ lhs } * rhs;
}
template <class T>
constexpr Rational<T> operator/(T lhs, const Rational<T>& rhs) noexcept {
	return Rational<T>{ lhs } / rhs;
}

namespace impl {
	template <class T>
	constexpr std::pair<T, T> CommonNumerators(const Rational<T>& lhs, const Rational<T>& rhs) noexcept {
		const auto common = std::lcm(lhs.Denominator(), rhs.Denominator());
		const auto lhsNum = lhs.Numerator() * (common / lhs.Denominator());
		const auto rhsNum = rhs.Numerator() * (common / rhs.Denominator());
		return { lhsNum, rhsNum };
	}
} // namespace impl

template <class T>
constexpr bool operator==(const Rational<T>& lhs, const Rational<T>& rhs) noexcept {
	const auto [lhsNum, rhsNum] = impl::CommonNumerators(lhs, rhs);
	return lhsNum == rhsNum;
}
template <class T>
constexpr bool operator!=(const Rational<T>& lhs, const Rational<T>& rhs) noexcept {
	return !(lhs == rhs);
}
template <class T>
constexpr bool operator<(const Rational<T>& lhs, const Rational<T>& rhs) noexcept {
	const auto [lhsNum, rhsNum] = impl::CommonNumerators(lhs, rhs);
	return lhsNum < rhsNum;
}
template <class T>
constexpr bool operator>(const Rational<T>& lhs, const Rational<T>& rhs) noexcept {
	const auto [lhsNum, rhsNum] = impl::CommonNumerators(lhs, rhs);
	return lhsNum > rhsNum;
}
template <class T>
constexpr bool operator<=(const Rational<T>& lhs, const Rational<T>& rhs) noexcept {
	return !(lhs > rhs);
}
template <class T>
constexpr bool operator>=(const Rational<T>& lhs, const Rational<T>& rhs) noexcept {
	return !(lhs < rhs);
}


} // namespace dspbb