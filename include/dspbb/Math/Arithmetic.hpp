#pragma once

#include <xsimd/xsimd.hpp>


namespace dspbb {


template <class T>
struct is_vectorized {
	static constexpr bool value = xsimd::simd_traits<T>::size > 1;
};


//------------------------------------------------------------------------------
// Generic calculation helpers.
//------------------------------------------------------------------------------

//--------------------------------------
// Vector-vector
//--------------------------------------
template <class R, class T, class U, class Op>
void Calculate(R* out, const T* a, const U* b, size_t length, Op op) {
	const R* last = out + length;
	for (; out != last; ++out, ++a, ++b) {
		*out = R(op(*a, *b));
	}
}


template <class T, class Op, std::enable_if_t<is_vectorized<T>::value, int> = 0>
void Calculate(T* out, const T* a, const T* b, size_t length, Op op) {
	using V = xsimd::simd_type<T>;
	constexpr size_t vsize = xsimd::simd_traits<T>::size;

	const size_t vlength = (length / vsize) * vsize;
	
	const T* vlast = out + vlength;
	for (; out < vlast; out += vsize, a += vsize, b += vsize) {
		V va, vb;
		va.load_unaligned(a);
		vb.load_unaligned(b);
		auto vr = op(va, vb);
		vr.store_unaligned(out);
	}

	Calculate<T, T, T>(out, a, b, length - vlength, op);
}


//--------------------------------------
// Scalar-vector
//--------------------------------------
template <class R, class T, class U, class Op, std::enable_if_t<!std::is_pointer<T>::value, int> = 0>
void Calculate(R* out, const T& a, const U* b, size_t length, Op op) {
	const R* last = out + length;
	for (; out != last; ++out, ++b) {
		*out = R(op(a, *b));
	}
}


template <class T, class Op, std::enable_if_t<!std::is_pointer<T>::value && is_vectorized<T>::value, int> = 0>
void Calculate(T* out, const T& a, const T* b, size_t length, Op op) {
	using V = xsimd::simd_type<T>;
	constexpr size_t vsize = xsimd::simd_traits<T>::size;

	const size_t vlength = (length / vsize) * vsize;

	const T* vlast = out + vlength;
	V va{ a };
	for (; out < vlast; out += vsize, b += vsize) {
		V vb;
		vb.load_unaligned(b);
		auto vr = op(va, vb);
		vr.store_unaligned(out);
	}

	Calculate<T, T, T>(out, a, b, length - vlength, op);
}


//--------------------------------------
// Vector-scalar
//--------------------------------------

template <class R, class T, class U, class Op, std::enable_if_t<!std::is_pointer<T>::value, int> = 0>
void Calculate(R* out, const T* a, const U& b, size_t length, Op op) {
	const R* last = out + length;
	for (; out != last; ++out, ++a) {
		*out = R(op(*a, b));
	}
}


template <class T, class Op, std::enable_if_t<!std::is_pointer<T>::value && is_vectorized<T>::value, int> = 0>
void Calculate(T* out, const T* a, const T& b, size_t length, Op op) {
	using V = xsimd::simd_type<T>;
	constexpr size_t vsize = xsimd::simd_traits<T>::size;

	const size_t vlength = (length / vsize) * vsize;

	const T* vlast = out + vlength;
	V vb{ b };
	for (; out < vlast; out += vsize, a += vsize) {
		V va;
		va.load_unaligned(a);
		auto vr = op(va, vb);
		vr.store_unaligned(out);
	}

	Calculate<T, T, T>(out, a, b, length - vlength, op);
}


//------------------------------------------------------------------------------
// Vector-vector operations.
//------------------------------------------------------------------------------

template <class R, class T, class U>
void Multiply(R* out, const T* a, const U* b, size_t length) {
	Calculate(out, a, b, length, [](const auto& a, const auto& b) { return a * b; });
}


template <class R, class T, class U>
void Divide(R* out, const T* a, const U* b, size_t length) {
	Calculate(out, a, b, length, [](const auto& a, const auto& b) { return a / b; });
}


template <class R, class T, class U>
void Add(R* out, const T* a, const U* b, size_t length) {
	Calculate(out, a, b, length, [](const auto& a, const auto& b) { return a + b; });
}


template <class R, class T, class U>
void Subtract(R* out, const T* a, const U* b, size_t length) {
	Calculate(out, a, b, length, [](const auto& a, const auto& b) { return a - b; });
}

//------------------------------------------------------------------------------
// Vector-scalar & scalar-vector operations.
//------------------------------------------------------------------------------

template <class R, class T, class U, std::enable_if_t<!std::is_pointer<U>::value, int> = 0>
void Multiply(R* out, const T* a, const U& b, size_t length) {
	Calculate(out, a, b, length, [](const auto& a, const auto& b) { return a * b; });
}


template <class R, class T, class U, std::enable_if_t<!std::is_pointer<U>::value, int> = 0>
void Divide(R* out, const T* a, const U& b, size_t length) {
	Calculate(out, a, b, length, [](const auto& a, const auto& b) { return a / b; });
}


template <class R, class T, class U, std::enable_if_t<!std::is_pointer<U>::value, int> = 0>
void Add(R* out, const T* a, const U& b, size_t length) {
	Calculate(out, a, b, length, [](const auto& a, const auto& b) { return a + b; });
}


template <class R, class T, class U, std::enable_if_t<!std::is_pointer<U>::value, int> = 0>
void Subtract(R* out, const T* a, const U& b, size_t length) {
	Calculate(out, a, b, length, [](const auto& a, const auto& b) { return a - b; });
}


template <class R, class T, class U, std::enable_if_t<!std::is_pointer<T>::value, int> = 0>
void Multiply(R* out, const T& a, const U* b, size_t length) {
	Calculate(out, a, b, length, [](const auto& a, const auto& b) { return a * b; });
}


template <class R, class T, class U, std::enable_if_t<!std::is_pointer<T>::value, int> = 0>
void Divide(R* out, const T& a, const U* b, size_t length) {
	Calculate(out, a, b, length, [](const auto& a, const auto& b) { return a / b; });
}


template <class R, class T, class U, std::enable_if_t<!std::is_pointer<T>::value, int> = 0>
void Add(R* out, const T& a, const U* b, size_t length) {
	Calculate(out, a, b, length, [](const auto& a, const auto& b) { return a + b; });
}


template <class R, class T, class U, std::enable_if_t<!std::is_pointer<T>::value, int> = 0>
void Subtract(R* out, const T& a, const U* b, size_t length) {
	Calculate(out, a, b, length, [](const auto& a, const auto& b) { return a - b; });
}

} // namespace dspbb