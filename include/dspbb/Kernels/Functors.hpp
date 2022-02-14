#pragma once

#include "../Utility/TypeTraits.hpp"

namespace dspbb {

//------------------------------------------------------------------------------
// Compensated operators
//------------------------------------------------------------------------------

struct compensated_operator_tag {};

template <class Operator>
struct is_operator_compensated : std::is_base_of<compensated_operator_tag, Operator> {};

template <class Operator>
constexpr bool is_operator_compensated_v = is_operator_compensated<Operator>::value;

template <class T, std::enable_if_t<std::is_convertible_v<int, T>, int> = 0>
T make_zero() {
	return static_cast<T>(0);
}

template <class T, class U = typename T::value_type, std::enable_if_t<!std::is_convertible_v<int, T>, int> = 0>
T make_zero() {
	return T{ static_cast<U>(0) };
}

template <class T = void>
struct plus_compensated : compensated_operator_tag {
	inline constexpr T operator()(const T& lhs, const T& rhs) const {
		return lhs + rhs;
	}
	inline constexpr T make_carry(const T&) const {
		return make_zero<T>();
	}
	inline constexpr T operator()(T& carry, const T& sum, const T& item) const {
		const T y = item - carry;
		const T t = sum + y;
		carry = (t - sum) - y;
		sum = t;
		return sum;
	}
};

template <>
struct plus_compensated<void> : compensated_operator_tag {
	template <class T, class U>
	inline constexpr auto operator()(T&& lhs, U&& rhs) const -> plus_result_t<T, U> {
		return std::forward<T>(lhs) + std::forward<U>(rhs);
	}
	template <class T, class U>
	inline constexpr auto make_carry(const plus_result_t<T, U>&) const -> plus_result_t<T, U> {
		return make_zero<T>();
	}
	template <class T, class U>
	inline constexpr auto operator()(plus_result_t<T, U>& carry, const T& sum, const U& item) const -> plus_result_t<T, U> {
		const auto y = item - carry;
		const auto t = sum + y;
		carry = (t - sum) - y;
		return t;
	}
};

//------------------------------------------------------------------------------
// Scalar-vector helpers
//------------------------------------------------------------------------------

template <class T>
struct multiplies_scalar_left {
	explicit multiplies_scalar_left(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<T, U> {
		return scalar * std::forward<U>(arg);
	}
	T scalar;
};

template <class T>
struct multiplies_scalar_right {
	explicit multiplies_scalar_right(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<U, T> {
		return std::forward<U>(arg) * scalar;
	}
	T scalar;
};

template <class T>
struct divides_scalar_left {
	explicit divides_scalar_left(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<T, U> {
		return scalar / std::forward<U>(arg);
	}
	T scalar;
};

template <class T>
struct divides_scalar_right {
	explicit divides_scalar_right(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<T, U> {
		return std::forward<U>(arg) / scalar;
	}
	T scalar;
};

template <class T>
struct plus_scalar_left {
	explicit plus_scalar_left(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<T, U> {
		return scalar + std::forward<U>(arg);
	}
	T scalar;
};

template <class T>
struct plus_scalar_right {
	explicit plus_scalar_right(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<U, T> {
		return std::forward<U>(arg) + scalar;
	}
	T scalar;
};

template <class T>
struct minus_scalar_left {
	explicit minus_scalar_left(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<T, U> {
		return scalar - std::forward<U>(arg);
	}
	T scalar;
};

template <class T>
struct minus_scalar_right {
	explicit minus_scalar_right(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<T, U> {
		return std::forward<U>(arg) - scalar;
	}
	T scalar;
};


} // namespace dspbb