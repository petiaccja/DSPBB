#pragma once

#pragma warning(push)
#pragma warning(disable : 4800)
#include <xsimd/xsimd.hpp>
#pragma warning(pop)

#include <numeric>


namespace dspbb::kernels {


template <class T>
struct is_vectorized {
	static constexpr bool value = xsimd::simd_traits<T>::size > 1;
};

template <class R, class T, class Op>
struct is_unary_vectorized {
	constexpr static bool get(...) { return false; }
	template <class R_ = R,
			  class T_ = T,
			  class Op_ = Op,
			  std::enable_if_t<(xsimd::simd_traits<T_>::size > 1)
								   && xsimd::simd_traits<R_>::size == xsimd::simd_traits<T_>::size
								   && std::is_same_v<xsimd::simd_type<R_>, std::invoke_result_t<Op, xsimd::simd_type<T_>>>,
							   int> = 0>
	constexpr static bool get(int) { return true; }
	static constexpr bool value = get(0);
};

template <class R, class T, class Op>
struct is_reduce_vectorized {
	constexpr static bool get(...) { return false; }
	template <class R_ = R,
			  class T_ = T,
			  class Op_ = Op,
			  std::enable_if_t<(xsimd::simd_traits<T_>::size > 1)
								   && xsimd::simd_traits<R_>::size == xsimd::simd_traits<T_>::size
								   && std::is_convertible_v<std::invoke_result_t<Op, xsimd::simd_type<R_>, xsimd::simd_type<T_>>, xsimd::simd_type<R_>>,
							   int> = 0>
	constexpr static bool get(int) { return true; }
	static constexpr bool value = get(0);
};

template <class R, class T, class ReduceOp, class MapOp>
struct is_map_reduce_vectorized {
	constexpr static bool get(...) { return false; }
	template <class R_ = R,
			  class T_ = T,
			  class ReduceOp_ = ReduceOp,
			  class MapOp_ = MapOp,
			  std::enable_if_t<(xsimd::simd_traits<T_>::size > 1)
								   && xsimd::simd_traits<R_>::size == xsimd::simd_traits<T_>::size
								   && std::is_convertible_v<std::invoke_result_t<ReduceOp, xsimd::simd_type<R_>, std::invoke_result_t<MapOp, xsimd::simd_type<T_>>>, xsimd::simd_type<R_>>,
							   int> = 0>
	constexpr static bool get(int) { return true; }
	static constexpr bool value = get(0);
};

template <class R, class T, class U, class ProductOp, class ReduceOp>
struct is_inner_product_vectorized {
	constexpr static bool get(...) { return false; }
	template <class R_ = R,
			  class T_ = T,
			  class U_ = U,
			  class ProductOp_ = ProductOp,
			  class ReduceOp_ = ReduceOp,
			  std::enable_if_t<(xsimd::simd_traits<T_>::size > 1)
								   && xsimd::simd_traits<R_>::size == xsimd::simd_traits<T_>::size
								   && xsimd::simd_traits<R_>::size == xsimd::simd_traits<U_>::size
								   && std::is_convertible_v<std::invoke_result_t<ReduceOp, xsimd::simd_type<R_>, std::invoke_result_t<ProductOp_, xsimd::simd_type<T_>, xsimd::simd_type<U_>>>, xsimd::simd_type<R_>>,
							   int> = 0>
	constexpr static bool get(int) { return true; }
	static constexpr bool value = get(0);
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

template <class R, class T, class U, class Op, std::enable_if_t<!std::is_pointer_v<T>, int> = 0>
void BinaryOperation(R* out, const T& a, const U* b, size_t length, Op op) {
	const R* last = out + length;
	for (; out != last; ++out, ++b) {
		*out = R(op(a, *b));
	}
}

template <class R, class T, class U, class Op, std::enable_if_t<!std::is_pointer_v<T>, int> = 0>
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

template <class T, class Op, std::enable_if_t<!std::is_pointer_v<T> && is_vectorized<T>::value, int> = 0>
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

template <class T, class Op, std::enable_if_t<!std::is_pointer_v<T> && is_vectorized<T>::value, int> = 0>
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

template <class R, class T, class U, class Op, std::enable_if_t<!std::is_pointer_v<T>, int> = 0>
void BinaryOperationVectorized(R* out, const T& a, const U* b, size_t length, Op op) {
	BinaryOperation(out, a, b, length, op);
}

template <class R, class T, class U, class Op, std::enable_if_t<!std::is_pointer_v<T>, int> = 0>
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

template <class R, class T, class Op, std::enable_if_t<!is_unary_vectorized<R, T, Op>::value, int> = 0>
void UnaryOperationVectorized(R* out, const T* in, size_t length, Op op) {
	return UnaryOperation(out, in, length, op);
}

template <class R, class T, class Op, std::enable_if_t<is_unary_vectorized<R, T, Op>::value, int> = 0>
void UnaryOperationVectorized(R* out, const T* in, size_t length, Op op) {
	using TV = xsimd::simd_type<T>;
	constexpr size_t vsize = xsimd::simd_traits<T>::size;

	const size_t vlength = (length / vsize) * vsize;

	const R* vlast = out + vlength;
	for (; out < vlast; out += vsize, in += vsize) {
		TV vin;
		vin.load_unaligned(in);
		const auto vr = op(vin);
		vr.store_unaligned(out);
	}

	UnaryOperation(out, in, length - vlength, op);
}

//------------------------------------------------------------------------------
// Reduction.
//------------------------------------------------------------------------------

template <class R, class T, class Op>
R Reduce(const T* in, size_t length, R init, Op op) {
	// TODO: use C++17 std::reduce
	return std::accumulate(in, in + length, init, op);
}

template <class R, class T, class Op, std::enable_if_t<!is_reduce_vectorized<R, T, Op>::value, int> = 0>
R ReduceVectorized(const T* in, size_t length, R init, Op op) {
	return Reduce(in, length, init, op);
}

template <class R, class T, class Op, std::enable_if_t<is_reduce_vectorized<R, T, Op>::value, int> = 0>
R ReduceVectorized(const T* in, size_t length, R init, Op op) {
	using TV = xsimd::simd_type<T>;
	using RV = xsimd::simd_type<R>;
	constexpr size_t vsize = xsimd::simd_traits<T>::size;

	const size_t vlength = (length / vsize) * vsize;
	const R* vlast = in + vlength;
	if (vlength > 0) {
		RV acc;
		acc.load_unaligned(in);
		in += vsize;
		for (; in < vlast; in += vsize) {
			TV vin;
			vin.load_unaligned(in);
			acc = op(acc, vin);
		}
		alignas(alignof(RV)) std::array<R, vsize> arr;
		acc.store_aligned(arr.data());
		init = Reduce(arr.data(), arr.size(), init, op);
	}

	return Reduce(in, length - vlength, init, op);
}

template <class R, class T, class ReduceOp, class MapOp>
R MapReduce(const T* in, size_t length, R init, ReduceOp reduceOp, MapOp mapOp) {
	// TODO: use C++17 std::transform_reduce
	R acc = std::move(init);
	const T* end = in + length;
	for (; in < end; ++in) {
		acc = reduceOp(acc, mapOp(*in));
	}
	return acc;
}

template <class R, class T, class ReduceOp, class MapOp, std::enable_if_t<!is_map_reduce_vectorized<R, T, ReduceOp, MapOp>::value, int> = 0>
R MapReduceVectorized(const T* in, size_t length, R init, ReduceOp reduceOp, MapOp mapOp) {
	return MapReduce(in, length, init, reduceOp, mapOp);
}

template <class R, class T, class ReduceOp, class MapOp, std::enable_if_t<is_map_reduce_vectorized<R, T, ReduceOp, MapOp>::value, int> = 0>
R MapReduceVectorized(const T* in, size_t length, R init, ReduceOp reduceOp, MapOp mapOp) {
	using TV = xsimd::simd_type<T>;
	using RV = xsimd::simd_type<R>;
	constexpr size_t vsize = xsimd::simd_traits<T>::size;

	const size_t vlength = (length / vsize) * vsize;
	const T* vlast = in + vlength;
	if (vlength > 0) {
		TV unmapped;
		unmapped.load_unaligned(in);
		RV acc = mapOp(unmapped);
		in += vsize;
		for (; in < vlast; in += vsize) {
			TV vin;
			vin.load_unaligned(in);
			acc = reduceOp(acc, mapOp(vin));
		}
		alignas(alignof(RV)) std::array<R, vsize> arr;
		acc.store_aligned(arr.data());
		init = Reduce(arr.data(), arr.size(), init, reduceOp);
	}

	return MapReduce(in, length - vlength, init, reduceOp, mapOp);
}


//------------------------------------------------------------------------------
// Inner product
//------------------------------------------------------------------------------

template <class R, class T, class U, class ProductOp, class ReduceOp>
R InnerProduct(const T* a, const U* b, size_t length, R init, ProductOp productOp, ReduceOp reduceOp) {
	R acc = std::move(init);
	for (size_t i = 0; i < length; ++i, ++a, ++b) {
		acc = reduceOp(acc, productOp(*a, *b));
	}
	return acc;
}

template <class R, class T, class U, class ProductOp, class ReduceOp, std::enable_if_t<!is_inner_product_vectorized<R, T, U, ProductOp, ReduceOp>::value, int> = 0>
R InnerProductVectorized(const T* a, const U* b, size_t length, R init, ProductOp productOp, ReduceOp reduceOp) {
	return InnerProduct(a, b, length, init, productOp, reduceOp);
}

template <class R, class T, class U, class ProductOp, class ReduceOp, std::enable_if_t<is_inner_product_vectorized<R, T, U, ProductOp, ReduceOp>::value, int> = 0>
R InnerProductVectorized(const T* a, const U* b, size_t length, R init, ProductOp productOp, ReduceOp reduceOp) {
	using TV = xsimd::simd_type<T>;
	using UV = xsimd::simd_type<T>;
	using RV = xsimd::simd_type<R>;
	constexpr size_t vsize = xsimd::simd_traits<T>::size;

	const size_t vlength = (length / vsize) * vsize;
	const T* alast = a + vlength;
	if (vlength > 0) {
		TV av;
		UV bv;
		av.load_unaligned(a);
		bv.load_unaligned(b);
		RV acc = productOp(av, bv);
		a += vsize;
		b += vsize;
		for (; a < alast; a += vsize, b += vsize) {
			av.load_unaligned(a);
			bv.load_unaligned(b);
			acc = reduceOp(acc, productOp(av, bv));
		}
		alignas(alignof(RV)) std::array<R, vsize> arr;
		acc.store_aligned(arr.data());
		init = Reduce(arr.data(), arr.size(), init, reduceOp);
	}

	return InnerProduct(a, b, length - vlength, init, productOp, reduceOp);
}



} // namespace dspbb
