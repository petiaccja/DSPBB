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

#include <numeric>


namespace dspbb::kernels {


// deprecated
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
								   && std::is_convertible_v<std::invoke_result_t<Op, xsimd::simd_type<T_>>, xsimd::simd_type<R_>>,
							   int> = 0>
	constexpr static bool get(int) { return true; }
	static constexpr bool value = get(0);
};



template <class T, class U, class UnaryOp>
struct is_transform_vectorized_1 {
	constexpr static bool get(...) { return false; }
	template <class T_ = T, class U_ = U, class UnaryOp_ = UnaryOp,
			  std::enable_if_t<(xsimd::simd_traits<T_>::size > 1)
								   && xsimd::simd_traits<T_>::size == xsimd::simd_traits<U_>::size
								   && std::is_convertible_v<std::invoke_result_t<UnaryOp, xsimd::simd_type<T_>>, xsimd::simd_type<U_>>,
							   int> = 0>
	constexpr static bool get(int) { return true; }
	static constexpr bool value = get(0);
};

template <class T1, class T2, class U, class BinaryOp>
struct is_transform_vectorized_2 {
	constexpr static bool get(...) { return false; }
	template <class T1_ = T1, class T2_ = T2, class U_ = U, class BinaryOp_ = BinaryOp,
			  std::enable_if_t<(xsimd::simd_traits<T1_>::size > 1)
								   && xsimd::simd_traits<T1_>::size == xsimd::simd_traits<U_>::size
								   && xsimd::simd_traits<T2_>::size == xsimd::simd_traits<U_>::size
								   && std::is_convertible_v<std::invoke_result_t<BinaryOp, xsimd::simd_type<T1_>, xsimd::simd_type<T2_>>, xsimd::simd_type<U_>>,
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
// Binary operations -- DEPRECATED BY TRANSFORM
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
// Utility
//------------------------------------------------------------------------------

template <class T, size_t N>
inline xsimd::batch<T, N> Load(const xsimd::batch<T, N>* p) {
	return xsimd::load_unaligned(reinterpret_cast<const T*>(p));
}

template <class T>
inline T Load(const T* p) {
	return *p;
}

//------------------------------------------------------------------------------
// Transform.
//------------------------------------------------------------------------------

template <class InputIter, class OutputIter, class UnaryOp>
inline auto Transform(InputIter first, InputIter last, OutputIter out, UnaryOp unaryOp)
	-> std::enable_if_t<is_random_access_iterator_v<InputIter> && is_random_access_iterator_v<OutputIter>, OutputIter> {
	using T = typename std::iterator_traits<InputIter>::value_type;
	using U = typename std::iterator_traits<OutputIter>::value_type;
	const auto count = std::distance(first, last);
	const T* pfirst = std::addressof(*first);
	const T* plast = pfirst + count;
	U* pout = std::addressof(*out);

	if constexpr (is_transform_vectorized_1<T, U, UnaryOp>::value) {
		constexpr size_t vectorWidth = xsimd::simd_traits<T>::size;

		const size_t vectorCount = count / vectorWidth;
		const auto* vectorLast = pfirst + vectorCount * vectorWidth;
		for (; pfirst != vectorLast; pfirst += vectorWidth, pout += vectorWidth) {
			xsimd::store_unaligned(pout, unaryOp(xsimd::load_unaligned(pfirst)));
		}
	}
	for (; pfirst != plast; ++pfirst, ++pout) {
		*pout = unaryOp(*pfirst);
	}
	return out + count;
}

template <class InputIter1, class InputIter2, class OutputIter, class BinaryOp>
inline auto Transform(InputIter1 first1, InputIter1 last1, InputIter2 first2, OutputIter out, BinaryOp binaryOp)
	-> std::enable_if_t<is_random_access_iterator_v<InputIter1> && is_random_access_iterator_v<InputIter2> && is_random_access_iterator_v<OutputIter>, OutputIter> {
	using T1 = typename std::iterator_traits<InputIter1>::value_type;
	using T2 = typename std::iterator_traits<InputIter2>::value_type;
	using U = typename std::iterator_traits<OutputIter>::value_type;
	const auto count = std::distance(first1, last1);
	const T1* pfirst1 = std::addressof(*first1);
	const T1* plast1 = pfirst1 + count;
	const T2* pfirst2 = std::addressof(*first2);
	U* pout = std::addressof(*out);

	if constexpr (is_transform_vectorized_2<T1, T2, U, BinaryOp>::value) {
		constexpr size_t vectorWidth = xsimd::simd_traits<T1>::size;

		const size_t vectorCount = count / vectorWidth;
		const auto* vectorLast = pfirst1 + vectorCount * vectorWidth;
		for (; pfirst1 != vectorLast; pfirst1 += vectorWidth, pfirst2 += vectorWidth, pout += vectorWidth) {
			xsimd::store_unaligned(pout, binaryOp(xsimd::load_unaligned(pfirst1), xsimd::load_unaligned(pfirst2)));
		}
	}
	for (; pfirst1 != plast1; ++pfirst1, ++pfirst2, ++pout) {
		*pout = binaryOp(*pfirst1, *pfirst2);
	}
	return out + count;
}


//------------------------------------------------------------------------------
// Reduce.
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

template <class T, class Init, class ReduceOp>
inline auto ReduceExplicit(const T* first, const T* last, const Init& init, ReduceOp reduceOp) -> Init {
	const size_t count = std::distance(first, last);
	const bool singlet = (count & 1) != 0;
	const bool doublet = (count & 2) != 0;
	const bool quadruplet = (count & 4) != 0;

	Init acc = init;
	if (singlet) {
		const auto val0 = Load(first);
		acc = reduceOp(acc, val0);
		first += 1;
	}
	if (doublet) {
		const auto val0 = Load(first);
		const auto val1 = Load(first + 1);
		acc = reduceOp(acc, reduceOp(val0, val1));
		first += 2;
	}
	if (quadruplet) {
		const auto val0 = Load(first);
		const auto val1 = Load(first + 1);
		const auto val2 = Load(first + 2);
		const auto val3 = Load(first + 3);
		acc = reduceOp(acc, reduceOp(reduceOp(val0, val1), reduceOp(val2, val3)));
		first += 4;
	}

	[[maybe_unused]] auto carry = make_compensation_carry<Init, T>(reduceOp, init);
	for (; first != last; first += 8) {
		const auto val0 = Load(first);
		const auto val1 = Load(first + 1);
		const auto val2 = Load(first + 2);
		const auto val3 = Load(first + 3);
		const auto val4 = Load(first + 4);
		const auto val5 = Load(first + 5);
		const auto val6 = Load(first + 6);
		const auto val7 = Load(first + 7);
		const auto partial = reduceOp(reduceOp(reduceOp(val0, val1), reduceOp(val2, val3)), reduceOp(reduceOp(val4, val5), reduceOp(val6, val7)));
		if constexpr (!is_operator_compensated_v<ReduceOp>) {
			acc = reduceOp(acc, partial);
		}
		else {
			acc = reduceOp(carry, acc, partial);
		}
	}
	return acc;
}

template <class Iter, class Init, class ReduceOp>
auto Reduce(Iter first, Iter last, Init init, ReduceOp reduceOp)
	-> std::enable_if_t<is_random_access_iterator_v<Iter>, Init> {
	using T = typename std::iterator_traits<Iter>::value_type;
	const auto count = std::distance(first, last);
	const T* pfirst = std::addressof(*first);
	const T* plast = pfirst + count;

	if constexpr (is_reduce_vectorized<Init, T, ReduceOp>::value) {
		constexpr size_t vectorWidth = xsimd::simd_traits<T>::size;

		const size_t vectorCount = count / vectorWidth;
		if (vectorCount != 0) {
			const auto vectorFirst = reinterpret_cast<const xsimd::simd_type<T>*>(pfirst);
			const auto vectorResult = ReduceExplicit(vectorFirst + 1, vectorFirst + vectorCount, Load(vectorFirst), reduceOp);
			pfirst += vectorCount * vectorWidth;
			init = ReduceBatch(vectorResult, std::move(init), reduceOp);
		}
		return std::reduce(pfirst, plast, init, reduceOp);
	}
	return ReduceExplicit(pfirst, plast, init, reduceOp);
}

//------------------------------------------------------------------------------
// Transform reduce.
//------------------------------------------------------------------------------


template <class T, class Init, class ReduceOp, class TransformOp>
inline auto TransformReduceExplicit(const T* first, const T* last, const Init& init, ReduceOp reduceOp, TransformOp transformOp) -> Init {
	const size_t count = std::distance(first, last);
	const bool singlet = (count & 1) != 0;
	const bool doublet = (count & 2) != 0;
	const bool quadruplet = (count & 4) != 0;

	Init acc = init;
	if (singlet) {
		const auto val0 = transformOp(Load(first));
		acc = reduceOp(acc, val0);
		first += 1;
	}
	if (doublet) {
		const auto val0 = transformOp(Load(first));
		const auto val1 = transformOp(Load(first + 1));
		acc = reduceOp(acc, reduceOp(val0, val1));
		first += 2;
	}
	if (quadruplet) {
		const auto val0 = transformOp(Load(first));
		const auto val1 = transformOp(Load(first + 1));
		const auto val2 = transformOp(Load(first + 2));
		const auto val3 = transformOp(Load(first + 3));
		acc = reduceOp(acc, reduceOp(reduceOp(val0, val1), reduceOp(val2, val3)));
		first += 4;
	}

	[[maybe_unused]] auto carry = make_compensation_carry<Init, T>(reduceOp, init);
	for (; first != last; first += 8) {
		const auto val0 = transformOp(Load(first));
		const auto val1 = transformOp(Load(first + 1));
		const auto val2 = transformOp(Load(first + 2));
		const auto val3 = transformOp(Load(first + 3));
		const auto val4 = transformOp(Load(first + 4));
		const auto val5 = transformOp(Load(first + 5));
		const auto val6 = transformOp(Load(first + 6));
		const auto val7 = transformOp(Load(first + 7));
		const auto partial = reduceOp(reduceOp(reduceOp(val0, val1), reduceOp(val2, val3)), reduceOp(reduceOp(val4, val5), reduceOp(val6, val7)));
		if constexpr (!is_operator_compensated_v<ReduceOp>) {
			acc = reduceOp(acc, partial);
		}
		else {
			acc = reduceOp(carry, acc, partial);
		}
	}
	return acc;
}

template <class Iter, class Init, class ReduceOp, class TransformOp>
auto TransformReduce(Iter first, Iter last, Init init, ReduceOp reduceOp, TransformOp transformOp)
	-> std::enable_if_t<is_random_access_iterator_v<Iter>, Init> {
	using T = typename std::iterator_traits<Iter>::value_type;
	const auto count = std::distance(first, last);
	const T* pfirst = std::addressof(*first);
	const T* plast = pfirst + count;

	if constexpr (is_map_reduce_vectorized<Init, T, ReduceOp, TransformOp>::value) {
		constexpr size_t vectorWidth = xsimd::simd_traits<T>::size;

		const size_t vectorCount = count / vectorWidth;
		if (vectorCount != 0) {
			const auto vectorFirst = reinterpret_cast<const xsimd::simd_type<T>*>(pfirst);
			const auto vectorResult = TransformReduceExplicit(vectorFirst + 1, vectorFirst + vectorCount, transformOp(Load(vectorFirst)), reduceOp, transformOp);
			pfirst += vectorCount * vectorWidth;
			init = ReduceBatch(vectorResult, std::move(init), reduceOp);
		}
		return std::transform_reduce(pfirst, plast, init, reduceOp, transformOp);
	}
	return TransformReduceExplicit(pfirst, plast, init, reduceOp, transformOp);
}

//------------------------------------------------------------------------------
// Inner product
//------------------------------------------------------------------------------


template <class T1, class T2, class Init, class ReduceOp, class ProductOp>
inline auto InnerProductExplicit(const T1* first1, const T1* last1, T2* first2, const Init& init, ReduceOp reduceOp, ProductOp productOp) -> Init {
	const size_t count = std::distance(first1, last1);
	const bool singlet = (count & 1) != 0;
	const bool doublet = (count & 2) != 0;
	const bool quadruplet = (count & 4) != 0;

	Init acc = init;
	if (singlet) {
		const auto val0 = productOp(Load(first1), Load(first2));
		acc = reduceOp(acc, val0);
		first1 += 1;
		first2 += 1;
	}
	if (doublet) {
		const auto val0 = productOp(Load(first1), Load(first2));
		const auto val1 = productOp(Load(first1 + 1), Load(first2 + 1));
		acc = reduceOp(acc, reduceOp(val0, val1));
		first1 += 2;
		first2 += 2;
	}
	if (quadruplet) {
		const auto val0 = productOp(Load(first1), Load(first2));
		const auto val1 = productOp(Load(first1 + 1), Load(first2 + 1));
		const auto val2 = productOp(Load(first1 + 2), Load(first2 + 2));
		const auto val3 = productOp(Load(first1 + 3), Load(first2 + 3));
		acc = reduceOp(acc, reduceOp(reduceOp(val0, val1), reduceOp(val2, val3)));
		first1 += 4;
		first2 += 4;
	}

	[[maybe_unused]] auto carry = make_compensation_carry<Init, std::invoke_result_t<ProductOp, T1, T2>>(reduceOp, init);
	for (; first1 != last1; first1 += 8, first2 += 8) {
		const auto val0 = productOp(Load(first1), Load(first2));
		const auto val1 = productOp(Load(first1 + 1), Load(first2 + 1));
		const auto val2 = productOp(Load(first1 + 2), Load(first2 + 2));
		const auto val3 = productOp(Load(first1 + 3), Load(first2 + 3));
		const auto val4 = productOp(Load(first1 + 4), Load(first2 + 4));
		const auto val5 = productOp(Load(first1 + 5), Load(first2 + 5));
		const auto val6 = productOp(Load(first1 + 6), Load(first2 + 6));
		const auto val7 = productOp(Load(first1 + 7), Load(first2 + 7));
		const auto partial = reduceOp(reduceOp(reduceOp(val0, val1), reduceOp(val2, val3)), reduceOp(reduceOp(val4, val5), reduceOp(val6, val7)));
		if constexpr (!is_operator_compensated_v<ReduceOp>) {
			acc = reduceOp(acc, partial);
		}
		else {
			acc = reduceOp(carry, acc, partial);
		}
	}
	return acc;
}

template <class Iter1, class Iter2, class Init, class ReduceOp, class ProductOp>
auto InnerProduct(Iter1 first1, Iter1 last1, Iter2 first2, Init init, ReduceOp reduceOp, ProductOp productOp)
	-> std::enable_if_t<is_random_access_iterator_v<Iter1> && is_random_access_iterator_v<Iter2>, Init> {
	using T1 = typename std::iterator_traits<Iter1>::value_type;
	using T2 = typename std::iterator_traits<Iter2>::value_type;

	const auto count = std::distance(first1, last1);
	const T1* pfirst1 = std::addressof(*first1);
	const T1* plast1 = pfirst1 + count;
	const T2* pfirst2 = std::addressof(*first2);

	if constexpr (is_inner_product_vectorized<Init, T1, T2, ProductOp, ReduceOp>::value) {
		constexpr size_t vectorWidth = xsimd::simd_traits<T1>::size;

		const size_t vectorCount = count / vectorWidth;
		if (vectorCount != 0) {
			const auto vectorFirst1 = reinterpret_cast<const xsimd::simd_type<T1>*>(pfirst1);
			const auto vectorFirst2 = reinterpret_cast<const xsimd::simd_type<T2>*>(pfirst2);
			const auto vectorInit = productOp(Load(vectorFirst1), Load(vectorFirst2));
			const auto vectorResult = InnerProductExplicit(vectorFirst1 + 1, vectorFirst1 + vectorCount, vectorFirst2 + 1, vectorInit, reduceOp, productOp);
			pfirst1 += vectorCount * vectorWidth;
			pfirst2 += vectorCount * vectorWidth;
			init = ReduceBatch(vectorResult, std::move(init), reduceOp);
		}
		return std::inner_product(pfirst1, plast1, pfirst2, init, reduceOp, productOp);
	}
	return InnerProductExplicit(pfirst1, plast1, pfirst2, init, reduceOp, productOp);
}

} // namespace dspbb::kernels
