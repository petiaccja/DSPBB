#pragma once

#ifdef _MSC_VER
	#pragma warning(push)
	#pragma warning(disable : 4800 4244)
#endif
#include <xsimd/xsimd.hpp>
#ifdef _MSC_VER
	#pragma warning(pop)
#endif

#include "dspbb/Utility/TypeTraits.hpp"

#include <cassert>
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
// Common operators.
//------------------------------------------------------------------------------

struct compensated_operator_tag {};

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
	inline constexpr auto operator()(T&& lhs, U&& rhs) const -> sum_type_t<T, U> {
		return lhs + rhs;
	}
	template <class T, class U>
	inline constexpr auto make_carry(const sum_type_t<T, U>& init) const -> sum_type_t<T, U> {
		return init - init; // Type may not be constructable from literal zero.
	}
	template <class T, class U>
	inline constexpr auto operator()(sum_type_t<T, U>& carry, T&& sum, U&& item) const -> sum_type_t<T, U> {
		const auto y = item - carry;
		const auto t = sum + y;
		carry = (t - sum) - y;
		return t;
	}
};

template <class Operator>
struct is_operator_compensated : std::is_base_of<compensated_operator_tag, Operator> {};

template <class Operator>
constexpr bool is_operator_compensated_v = is_operator_compensated<Operator>::value;

template <class Arg1, class Arg2, class CarryT, class Operator>
inline auto make_compensation_carry(const Operator& op, const CarryT& init) -> std::invoke_result_t<decltype(&Operator::template make_carry<Arg1, Arg2>), Operator*, CarryT> {
	return op.template make_carry<Arg1, Arg2>(init);
}

template <class Arg1, class Arg2, class CarryT, class Operator>
inline auto make_compensation_carry(const Operator& op, const CarryT& init) -> std::invoke_result_t<decltype(&Operator::make_carry), Operator*, CarryT> {
	return op.make_carry(init);
}

template <class Arg1, class Arg2, class CarryT, class Operator>
inline auto make_compensation_carry(const Operator&, const CarryT&) -> std::enable_if_t<!is_operator_compensated_v<Operator>, compensated_operator_tag> {
	return compensated_operator_tag{}; // Return a useless tag just for the sake of compiling in generic contexts.
}


//------------------------------------------------------------------------------
// Reduction.
//------------------------------------------------------------------------------

template <class T, size_t N, class Init, class ReduceOp>
inline Init ReduceBatch(const xsimd::batch<T, N>& batch, Init init, ReduceOp reduceOp) {
	alignas(alignof(xsimd::batch<T, N>)) std::array<T, N> elements;
	batch.store_unaligned(elements.data());
	return std::reduce(elements.begin(), elements.end(), std::move(init), std::move(reduceOp));
}

template <class T, class U, size_t N, class Init>
inline Init ReduceBatch(const xsimd::batch<T, N>& batch, Init init, plus_compensated<U>) {
	return init + xsimd::hadd(batch);
}

template <class T, class U, size_t N, class Init>
inline Init ReduceBatch(const xsimd::batch<T, N>& batch, Init init, std::plus<U>) {
	return init + xsimd::hadd(batch);
}

template <class Init, class T, class ReduceOp>
inline auto ReduceExplicit(const T* first, const T* last, const Init& init, ReduceOp reduceOp) -> Init {
	const size_t count = std::distance(first, last);
	const bool singlet = (count & 1) != 0;
	const bool doublet = (count & 2) != 0;
	const bool quadruplet = (count & 4) != 0;

	Init acc = init;
	if (singlet) {
		acc = reduceOp(acc, first[0]);
		first += 1;
	}
	if (doublet) {
		acc = reduceOp(acc, reduceOp(first[0], first[1]));
		first += 2;
	}
	if (quadruplet) {
		acc = reduceOp(acc, reduceOp(reduceOp(first[0], first[1]), reduceOp(first[2], first[3])));
		first += 4;
	}

	[[maybe_unused]] auto carry = make_compensation_carry<Init, T>(reduceOp, init);
	for (; first != last; first += 8) {
		if constexpr (!is_operator_compensated_v<ReduceOp>) {
			acc = reduceOp(acc,
						   reduceOp(reduceOp(reduceOp(first[0], first[1]), reduceOp(first[2], first[3])),
									reduceOp(reduceOp(first[4], first[5]), reduceOp(first[6], first[7]))));
		}
		else {
			acc = reduceOp(carry, acc,
						   reduceOp(reduceOp(reduceOp(first[0], first[1]), reduceOp(first[2], first[3])),
									reduceOp(reduceOp(first[4], first[5]), reduceOp(first[6], first[7]))));
		}
	}
	return acc;
}

template <class Init, class T, class ReduceOp>
auto Reduce(const T* data, size_t count, Init init, ReduceOp reduceOp) -> Init {
	if constexpr (is_reduce_vectorized<Init, T, ReduceOp>::value) {
		constexpr size_t vectorWidth = xsimd::simd_traits<T>::size;

		const size_t vectorCount = count / vectorWidth;
		if (vectorCount != 0) {
			const auto vectorData = reinterpret_cast<const xsimd::simd_type<T>*>(data);
			const auto vectorResult = ReduceExplicit(vectorData + 1, vectorData + vectorCount, vectorData[0], reduceOp);
			data += vectorCount * vectorWidth;
			count -= vectorCount * vectorWidth;
			init = ReduceBatch(vectorResult, std::move(init), reduceOp);
		}
		return ReduceExplicit(data, data + count, init, reduceOp);
	}
	return ReduceExplicit(data, data + count, init, reduceOp);
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
	using UV = xsimd::simd_type<U>;
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



} // namespace dspbb::kernels
