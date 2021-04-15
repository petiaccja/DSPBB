#pragma once

#pragma warning(push)
#pragma warning(disable : 4800)
#include <xsimd/xsimd.hpp>
#pragma warning(pop)


namespace dspbb {


template <class T>
struct is_vectorized {
	static constexpr bool value = xsimd::simd_traits<T>::size > 1;
};


//------------------------------------------------------------------------------
// Binary operations.
//------------------------------------------------------------------------------

template <class R, class T, class U, class Op>
void BinaryOperation(R* out, const T* a, const U* b, size_t length, Op op) {
	const R* last = out + length;
	for (; out != last; ++out, ++a, ++b) {
		*out = R(op(*a, *b));
	}
}

template <class R, class T, class U, class Op, std::enable_if_t<!std::is_pointer<T>::value, int> = 0>
void BinaryOperation(R* out, const T& a, const U* b, size_t length, Op op) {
	const R* last = out + length;
	for (; out != last; ++out, ++b) {
		*out = R(op(a, *b));
	}
}

template <class R, class T, class U, class Op, std::enable_if_t<!std::is_pointer<T>::value, int> = 0>
void BinaryOperation(R* out, const T* a, const U& b, size_t length, Op op) {
	const R* last = out + length;
	for (; out != last; ++out, ++a) {
		*out = R(op(*a, b));
	}
}

template <class T, class Op, std::enable_if_t<is_vectorized<T>::value, int> = 0>
void BinaryOperationVectorized(T* out, const T* a, const T* b, size_t length, Op op) {
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

	BinaryOperation(out, a, b, length - vlength, op);
}

template <class T, class Op, std::enable_if_t<!std::is_pointer<T>::value && is_vectorized<T>::value, int> = 0>
void BinaryOperationVectorized(T* out, const T& a, const T* b, size_t length, Op op) {
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

	BinaryOperation(out, a, b, length - vlength, op);
}

template <class T, class Op, std::enable_if_t<!std::is_pointer<T>::value && is_vectorized<T>::value, int> = 0>
void BinaryOperationVectorized(T* out, const T* a, const T& b, size_t length, Op op) {
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

	BinaryOperation(out, a, b, length - vlength, op);
}


template <class R, class T, class U, class Op>
void BinaryOperationVectorized(R* out, const T* a, const U* b, size_t length, Op op) {
	BinaryOperation(out, a, b, length, op);
}

template <class R, class T, class U, class Op, std::enable_if_t<!std::is_pointer<T>::value, int> = 0>
void BinaryOperationVectorized(R* out, const T& a, const U* b, size_t length, Op op) {
	BinaryOperation(out, a, b, length, op);
}

template <class R, class T, class U, class Op, std::enable_if_t<!std::is_pointer<T>::value, int> = 0>
void BinaryOperationVectorized(R* out, const T* a, const U& b, size_t length, Op op) {
	BinaryOperation(out, a, b, length, op);
}

//------------------------------------------------------------------------------
// Unary operations.
//------------------------------------------------------------------------------

template <class R, class T, class Op>
void UnaryOperation(R* out, const T* in, size_t length, Op op) {
	const R* last = out + length;
	for (; out != last; ++out, ++in) {
		*out = R(op(*in));
	}
}

template <class R, class T, class Op>
void UnaryOperationVectorized(R* out, const T* in, size_t length, Op op) {
	return UnaryOperation(out, in, length, op);
}

template <class T, class Op, std::enable_if_t<!is_vectorized<T>::value, int> = 0>
void UnaryOperationVectorized(T* out, T* in, size_t length, Op op) {
	using V = xsimd::simd_type<T>;
	constexpr size_t vsize = xsimd::simd_traits<T>::size;

	const size_t vlength = (length / vsize) * vsize;

	const T* vlast = out + vlength;
	for (; out < vlast; out += vsize, in += vsize) {
		V vin;
		vin.load_unaligned(in);
		auto vr = op(vin);
		vr.store_unaligned(out);
	}

	UnaryOperation(out, in, length - vlength, op);
}

template <class R, class T, class Op, class VecOp>
void UnaryOperationVectorized(R* out, T* in, size_t length, size_t stride, VecOp vop, Op op) {
	const size_t vlength = (length / stride) * stride;

	const R* vlast = out + vlength;
	const R* last = out + length;

	for (; out < vlast; out += stride, in += stride) {
		vop(out, in);
	}
	for (; out < last; out += 1, in += 1) {
		op(out, in);
	}	
}


} // namespace dspbb
