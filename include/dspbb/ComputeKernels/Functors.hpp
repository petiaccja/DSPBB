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

template <class T = void>
struct plus_compensated : compensated_operator_tag {
	inline constexpr T operator()(const T& lhs, const T& rhs) const {
		return lhs + rhs;
	}
	inline constexpr T make_carry(const T& init) const {
		return init - init; // Type may not be constructable from literal zero.
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
		return lhs + rhs;
	}
	template <class T, class U>
	inline constexpr auto make_carry(const plus_result_t<T, U>& init) const -> plus_result_t<T, U> {
		return init - init; // Type may not be constructable from literal zero.
	}
	template <class T, class U>
	inline constexpr auto operator()(plus_result_t<T, U>& carry, T&& sum, U&& item) const -> plus_result_t<T, U> {
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
	multiplies_scalar_left(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<T, U> {
		return scalar * arg;
	}
	T scalar;
};

template <class T>
struct multiplies_scalar_right {
	multiplies_scalar_right(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<U, T> {
		return arg * scalar;
	}
	T scalar;
};

template <class T>
struct divides_scalar_left {
	divides_scalar_left(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<T, U> {
		return scalar / arg;
	}
	T scalar;
};

template <class T>
struct divides_scalar_right {
	divides_scalar_right(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<T, U> {
		return arg / scalar;
	}
	T scalar;
};

template <class T>
struct plus_scalar_left {
	plus_scalar_left(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<T, U> {
		return scalar + arg;
	}
	T scalar;
};

template <class T>
struct plus_scalar_right {
	plus_scalar_right(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<U, T> {
		return arg + scalar;
	}
	T scalar;
};

template <class T>
struct minus_scalar_left {
	minus_scalar_left(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<T, U> {
		return scalar - arg;
	}
	T scalar;
};

template <class T>
struct minus_scalar_right {
	minus_scalar_right(T scalar) : scalar(std::move(scalar)) {}
	template <class U>
	constexpr auto operator()(U&& arg) -> multiplies_result_t<T, U> {
		return arg - scalar;
	}
	T scalar;
};


} // namespace dspbb