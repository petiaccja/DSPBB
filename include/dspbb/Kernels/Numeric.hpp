#pragma once

#ifdef _MSC_VER
	#pragma warning(push)
	#pragma warning(disable : 4800 4244)
#endif
#include <xsimd/xsimd.hpp>
#ifdef _MSC_VER
	#pragma warning(pop)
#endif

#include "../Utility/TypeTraits.hpp"
#include "Functors.hpp"
#include "Utility.hpp"

#include <array>
#include <numeric>


namespace dspbb::kernels {

//------------------------------------------------------------------------------
// Vectorization possibility
//------------------------------------------------------------------------------

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
// Utility
//------------------------------------------------------------------------------

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
		using V = xsimd::batch<T>;
		using VU = xsimd::batch<U>;
		constexpr size_t vectorWidth = xsimd::simd_traits<T>::size;

		const size_t vectorCount = count / vectorWidth;
		const auto* vectorLast = pfirst + vectorCount * vectorWidth;
		for (; pfirst != vectorLast; pfirst += vectorWidth, pout += vectorWidth) {
			const VU result = unaryOp(V::load_unaligned(pfirst));
			result.store_unaligned(pout);
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
		using V1 = xsimd::batch<T1>;
		using V2 = xsimd::batch<T2>;
		using VU = xsimd::batch<U>;
		constexpr size_t vectorWidth = xsimd::simd_traits<T1>::size;

		const size_t vectorCount = count / vectorWidth;
		const auto* vectorLast = pfirst1 + vectorCount * vectorWidth;
		for (; pfirst1 != vectorLast; pfirst1 += vectorWidth, pfirst2 += vectorWidth, pout += vectorWidth) {
			const VU result = binaryOp(V1::load_unaligned(pfirst1), V2::load_unaligned(pfirst2));
			result.store_unaligned(pout);
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

template <class T, class Arch, class Init, class ReduceOp>
inline Init ReduceBatch(const xsimd::batch<T, Arch>& batch, Init init, ReduceOp reduceOp) {
	constexpr size_t batchSize = xsimd::revert_simd_traits<xsimd::batch<T, Arch>>::size;
	alignas(alignof(xsimd::batch<T, Arch>)) std::array<T, batchSize> elements;
	batch.store_unaligned(elements.data());
	return std::reduce(elements.begin(), elements.end(), std::move(init), std::move(reduceOp));
}

template <class T, class U, class Arch, class Init>
inline Init ReduceBatch(const xsimd::batch<T, Arch>& batch, Init init, plus_compensated<U>) {
	return init + xsimd::hadd(batch);
}

template <class T, class U, class Arch, class Init>
inline Init ReduceBatch(const xsimd::batch<T, Arch>& batch, Init init, std::plus<U>) {
	return init + xsimd::hadd(batch);
}

template <class T, class Init, class ReduceOp>
inline auto ReduceExplicit(const T* first, const T* last, const Init& init, ReduceOp reduceOp) -> Init {
	using V = std::conditional_t<xsimd::is_batch<Init>::value, xsimd::simd_type<T>, T>;
	constexpr size_t stride = xsimd::is_batch<Init>::value ? xsimd::revert_simd_traits<Init>::size : 1;
	const size_t count = std::distance(first, last) / stride;
	const bool singlet = (count & 1) != 0;
	const bool doublet = (count & 2) != 0;
	const bool quadruplet = (count & 4) != 0;

	Init acc = init;
	if (singlet) {
		const auto val0 = uniform_load_unaligned<V>(first);
		acc = reduceOp(acc, val0);
		first += 1 * stride;
	}
	if (doublet) {
		const auto val0 = uniform_load_unaligned<V>(first);
		const auto val1 = uniform_load_unaligned<V>(first + 1 * stride);
		acc = reduceOp(acc, reduceOp(val0, val1));
		first += 2 * stride;
	}
	if (quadruplet) {
		const auto val0 = uniform_load_unaligned<V>(first);
		const auto val1 = uniform_load_unaligned<V>(first + 1 * stride);
		const auto val2 = uniform_load_unaligned<V>(first + 2 * stride);
		const auto val3 = uniform_load_unaligned<V>(first + 3 * stride);
		acc = reduceOp(acc, reduceOp(reduceOp(val0, val1), reduceOp(val2, val3)));
		first += 4 * stride;
	}

	[[maybe_unused]] auto carry = make_compensation_carry<Init, T>(reduceOp, init);
	for (; first != last; first += 8 * stride) {
		const auto val0 = uniform_load_unaligned<V>(first);
		const auto val1 = uniform_load_unaligned<V>(first + 1 * stride);
		const auto val2 = uniform_load_unaligned<V>(first + 2 * stride);
		const auto val3 = uniform_load_unaligned<V>(first + 3 * stride);
		const auto val4 = uniform_load_unaligned<V>(first + 4 * stride);
		const auto val5 = uniform_load_unaligned<V>(first + 5 * stride);
		const auto val6 = uniform_load_unaligned<V>(first + 6 * stride);
		const auto val7 = uniform_load_unaligned<V>(first + 7 * stride);
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
		using V = const xsimd::simd_type<T>;
		constexpr size_t vectorWidth = xsimd::simd_traits<T>::size;

		const size_t vectorCount = count / vectorWidth;
		if (vectorCount != 0) {
			const auto vinit = uniform_load_unaligned<V>(pfirst);
			const auto vectorResult = ReduceExplicit(pfirst + vectorWidth, pfirst + vectorCount * vectorWidth, vinit, reduceOp);
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
	using V = std::conditional_t<xsimd::is_batch<Init>::value, xsimd::simd_type<T>, T>;
	constexpr size_t stride = xsimd::is_batch<Init>::value ? xsimd::revert_simd_traits<Init>::size : 1;
	const size_t count = std::distance(first, last) / stride;
	const bool singlet = (count & 1) != 0;
	const bool doublet = (count & 2) != 0;
	const bool quadruplet = (count & 4) != 0;

	Init acc = init;
	if (singlet) {
		const auto val0 = transformOp(uniform_load_unaligned<V>(first));
		acc = reduceOp(acc, val0);
		first += 1 * stride;
	}
	if (doublet) {
		const auto val0 = transformOp(uniform_load_unaligned<V>(first));
		const auto val1 = transformOp(uniform_load_unaligned<V>(first + 1 * stride));
		acc = reduceOp(acc, reduceOp(val0, val1));
		first += 2 * stride;
	}
	if (quadruplet) {
		const auto val0 = transformOp(uniform_load_unaligned<V>(first));
		const auto val1 = transformOp(uniform_load_unaligned<V>(first + 1 * stride));
		const auto val2 = transformOp(uniform_load_unaligned<V>(first + 2 * stride));
		const auto val3 = transformOp(uniform_load_unaligned<V>(first + 3 * stride));
		acc = reduceOp(acc, reduceOp(reduceOp(val0, val1), reduceOp(val2, val3)));
		first += 4 * stride;
	}

	[[maybe_unused]] auto carry = make_compensation_carry<Init, T>(reduceOp, init);
	for (; first != last; first += 8 * stride) {
		const auto val0 = transformOp(uniform_load_unaligned<V>(first));
		const auto val1 = transformOp(uniform_load_unaligned<V>(first + 1 * stride));
		const auto val2 = transformOp(uniform_load_unaligned<V>(first + 2 * stride));
		const auto val3 = transformOp(uniform_load_unaligned<V>(first + 3 * stride));
		const auto val4 = transformOp(uniform_load_unaligned<V>(first + 4 * stride));
		const auto val5 = transformOp(uniform_load_unaligned<V>(first + 5 * stride));
		const auto val6 = transformOp(uniform_load_unaligned<V>(first + 6 * stride));
		const auto val7 = transformOp(uniform_load_unaligned<V>(first + 7 * stride));
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
		using V = xsimd::simd_type<T>;
		constexpr size_t vectorWidth = xsimd::simd_traits<T>::size;

		const size_t vectorCount = count / vectorWidth;
		if (vectorCount != 0) {
			const auto vectorResult = TransformReduceExplicit(pfirst + vectorWidth, pfirst + vectorCount * vectorWidth, transformOp(uniform_load_unaligned<V>(pfirst)), reduceOp, transformOp);
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
inline auto InnerProductExplicit(const T1* first1, const T1* last1, const T2* first2, const Init& init, ReduceOp reduceOp, ProductOp productOp) -> Init {
	using V1 = std::conditional_t<xsimd::is_batch<Init>::value, xsimd::simd_type<T1>, T1>;
	using V2 = std::conditional_t<xsimd::is_batch<Init>::value, xsimd::simd_type<T2>, T2>;
	constexpr size_t stride = xsimd::is_batch<Init>::value ? xsimd::revert_simd_traits<Init>::size : 1;

	const size_t count = std::distance(first1, last1) / stride;
	const bool singlet = (count & 1) != 0;
	const bool doublet = (count & 2) != 0;
	const bool quadruplet = (count & 4) != 0;

	Init acc = init;
	if (singlet) {
		const auto val0 = productOp(uniform_load_unaligned<V1>(first1), uniform_load_unaligned<V2>(first2));
		acc = reduceOp(acc, val0);
		first1 += 1 * stride;
		first2 += 1 * stride;
	}
	if (doublet) {
		const auto val0 = productOp(uniform_load_unaligned<V1>(first1), uniform_load_unaligned<V2>(first2));
		const auto val1 = productOp(uniform_load_unaligned<V1>(first1 + 1 * stride), uniform_load_unaligned<V2>(first2 + 1 * stride));
		acc = reduceOp(acc, reduceOp(val0, val1));
		first1 += 2 * stride;
		first2 += 2 * stride;
	}
	if (quadruplet) {
		const auto val0 = productOp(uniform_load_unaligned<V1>(first1), uniform_load_unaligned<V2>(first2));
		const auto val1 = productOp(uniform_load_unaligned<V1>(first1 + 1 * stride), uniform_load_unaligned<V2>(first2 + 1 * stride));
		const auto val2 = productOp(uniform_load_unaligned<V1>(first1 + 2 * stride), uniform_load_unaligned<V2>(first2 + 2 * stride));
		const auto val3 = productOp(uniform_load_unaligned<V1>(first1 + 3 * stride), uniform_load_unaligned<V2>(first2 + 3 * stride));
		acc = reduceOp(acc, reduceOp(reduceOp(val0, val1), reduceOp(val2, val3)));
		first1 += 4 * stride;
		first2 += 4 * stride;
	}

	[[maybe_unused]] auto carry = make_compensation_carry<Init, std::invoke_result_t<ProductOp, V1, V2>>(reduceOp, init);
	for (; first1 != last1; first1 += 8 * stride, first2 += 8 * stride) {
		const auto val0 = productOp(uniform_load_unaligned<V1>(first1), uniform_load_unaligned<V2>(first2));
		const auto val1 = productOp(uniform_load_unaligned<V1>(first1 + 1 * stride), uniform_load_unaligned<V2>(first2 + 1 * stride));
		const auto val2 = productOp(uniform_load_unaligned<V1>(first1 + 2 * stride), uniform_load_unaligned<V2>(first2 + 2 * stride));
		const auto val3 = productOp(uniform_load_unaligned<V1>(first1 + 3 * stride), uniform_load_unaligned<V2>(first2 + 3 * stride));
		const auto val4 = productOp(uniform_load_unaligned<V1>(first1 + 4 * stride), uniform_load_unaligned<V2>(first2 + 4 * stride));
		const auto val5 = productOp(uniform_load_unaligned<V1>(first1 + 5 * stride), uniform_load_unaligned<V2>(first2 + 5 * stride));
		const auto val6 = productOp(uniform_load_unaligned<V1>(first1 + 6 * stride), uniform_load_unaligned<V2>(first2 + 6 * stride));
		const auto val7 = productOp(uniform_load_unaligned<V1>(first1 + 7 * stride), uniform_load_unaligned<V2>(first2 + 7 * stride));
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
		using V1 = xsimd::simd_type<T1>;
		using V2 = xsimd::simd_type<T2>;
		constexpr size_t vectorWidth = xsimd::simd_traits<T1>::size;

		const size_t vectorCount = count / vectorWidth;
		if (vectorCount != 0) {
			const auto vectorInit = productOp(uniform_load_unaligned<V1>(pfirst1), uniform_load_unaligned<V2>(pfirst2));
			const auto vectorResult = InnerProductExplicit(pfirst1 + 1 * vectorWidth, pfirst1 + vectorCount * vectorWidth, pfirst2 + vectorWidth, vectorInit, reduceOp, productOp);
			pfirst1 += vectorCount * vectorWidth;
			pfirst2 += vectorCount * vectorWidth;
			init = ReduceBatch(vectorResult, std::move(init), reduceOp);
		}
		return std::inner_product(pfirst1, plast1, pfirst2, init, reduceOp, productOp);
	}
	return InnerProductExplicit(pfirst1, plast1, pfirst2, init, reduceOp, productOp);
}

} // namespace dspbb::kernels
