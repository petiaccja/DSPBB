#pragma once

#include "../Utility/Interval.hpp"
#include "Math.hpp"
#include "Numeric.hpp"
#include "Utility.hpp"

#include <cassert>
#include <iterator>
#include <utility>


namespace dspbb::kernels {


template <class T1, class T2, class OutT>
struct is_convolution_reduce_vectorized {
	constexpr static bool get(...) { return false; }
	template <class T1_ = T1, class T2_ = T2, class OutT_ = OutT,
			  std::enable_if_t<(xsimd::simd_traits<OutT_>::size > 1)
								   && std::is_invocable_v<std::plus<>, xsimd::simd_type<OutT_>, std::invoke_result_t<std::multiplies<>, xsimd::simd_type<T1_>, xsimd::simd_type<T2_>>>,
							   int> = 0>
	constexpr static bool get(int) { return true; }
	static constexpr bool value = get(0);
};

template <class T>
struct convolution_fma {
	explicit convolution_fma(T multiplier) : multiplier(std::move(multiplier)) {}

	template <class U1, class U2, std::enable_if_t<!xsimd::is_batch<std::decay_t<U1>>::value || !xsimd::is_batch<std::decay_t<U2>>::value, int> = 0>
	constexpr auto operator()(U1&& accumulator, U2&& increase) const
		-> decltype(math_functions::fma(std::declval<U2>(), std::declval<T>(), std::declval<U1>())) {
		return math_functions::fma(increase, multiplier, accumulator);
	}

	template <class U1, class U2, std::enable_if_t<xsimd::is_batch<std::decay_t<U1>>::value && xsimd::is_batch<std::decay_t<U2>>::value, int> = 0>
	constexpr auto operator()(const U1& accumulator, const U2& increase) const
		-> decltype(math_functions::fma(std::declval<U2>(), std::declval<xsimd::simd_type<T>>(), std::declval<U1>())) {
		return math_functions::fma(increase, xsimd::simd_type<T>(multiplier), accumulator);
	}

	T multiplier;
};


template <class Iter1, class Iter2, class IterOut>
void ConvolutionNaive(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, IterOut firstOut, IterOut lastOut, ptrdiff_t n, bool accumulate = false) {
	const ptrdiff_t len1 = std::distance(first1, last1);
	const ptrdiff_t len2 = std::distance(first2, last2);

	for (; firstOut != lastOut; ++firstOut, ++n) {
		auto accumulator = accumulate ? *firstOut : typename std::iterator_traits<IterOut>::value_type(0);
		ptrdiff_t mLast = std::min(len1, n + 1);
		ptrdiff_t mFirst = std::max(ptrdiff_t(0), n - len2 + 1);
		for (; mFirst < mLast; ++mFirst) {
			accumulator += first1[mFirst] * first2[n - mFirst];
		}
		*firstOut = accumulator;
	}
}

template <class Iter1, class Iter2, class IterOut>
void ConvolutionSlide(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, IterOut firstOut, IterOut lastOut, ptrdiff_t n, bool accumulate = false) {
	const ptrdiff_t len1 = std::distance(first1, last1);
	const ptrdiff_t len2 = std::distance(first2, last2);
	const ptrdiff_t lenOut = std::distance(firstOut, lastOut);

	// We want input #2 to be at least say 512 for vectorization, but not more to keep it in L1 cache.
	if (std::min(len1, len2) > 512) {
		if (len1 < len2) {
			return ConvolutionSlide(first2, last2, first1, last1, firstOut, lastOut, n, accumulate);
		}
	}
	else {
		if (len2 < len1) {
			return ConvolutionSlide(first2, last2, first1, last1, firstOut, lastOut, n, accumulate);
		}
	}

	if (!accumulate) {
		std::fill(firstOut, lastOut, typename std::iterator_traits<IterOut>::value_type(0));
	}

	const Interval multiplierRange = { std::max(ptrdiff_t(0), n - len2 + 1), std::min(len1, n + lenOut) };
	const Interval outRange = { n, n + lenOut };
	Interval slidingRange = { multiplierRange.first, multiplierRange.first + len2 };
	for (; slidingRange.first < multiplierRange.last; slidingRange += ptrdiff_t(1)) {
		const auto multiplier = *(first1 + slidingRange.first);
		const auto writeRange = Intersection(outRange, slidingRange);
		const auto writeRangeOut = writeRange - n;
		const auto writeFirst2 = writeRange.first - slidingRange.first;

		Transform(firstOut + writeRangeOut.first, firstOut + writeRangeOut.last,
				  first2 + writeFirst2,
				  firstOut + writeRangeOut.first,
				  convolution_fma{ multiplier });
	}
}

template <class Iter1, class Iter2, class IterOut>
void ConvolutionReduce(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, IterOut firstOut, IterOut lastOut, ptrdiff_t n, bool accumulate = false) {
	assert(accumulate == false); // Accumulate is not implemented in this function yet, so better make sure it's not misused.

	const ptrdiff_t len1 = std::distance(first1, last1);
	const ptrdiff_t len2 = std::distance(first2, last2);

	constexpr ptrdiff_t vectorWidth = 8;
	using OutT = typename std::iterator_traits<IterOut>::value_type;

	for (; firstOut < lastOut; n += vectorWidth) {
		std::array<OutT, vectorWidth> accumulator;
		std::array<OutT, vectorWidth> data;
		std::fill(accumulator.begin(), accumulator.end(), OutT(0));
		std::fill(data.begin(), data.end(), OutT(0));

		ptrdiff_t mFirst = std::max(ptrdiff_t(0), n - len2 + 1);
		ptrdiff_t mLast = std::min(len1, n + vectorWidth);
		for (; mFirst < mLast; ++mFirst) {
			for (ptrdiff_t i = 0; i < vectorWidth; ++i) {
				const auto access = n + i - mFirst;
				data[i] = access >= 0 && access < len2 ? first2[access] : 0;
			}

			for (ptrdiff_t i = 0; i < vectorWidth; ++i) {
				accumulator[i] += first1[mFirst] * data[i];
			}
		}

		for (size_t i = 0; i < vectorWidth && firstOut < lastOut; ++i, ++firstOut) {
			*firstOut = accumulator[i];
		}
	}
}


template <bool Vectorize, class Iter1, class Iter2, class OutV, class ReduceOp>
OutV ConvolutionReduceLoop(Iter1 first1, Iter2 first2, OutV init, ptrdiff_t count, ReduceOp reduceOp) {
	using T1 = typename std::iterator_traits<Iter1>::value_type;
	using T2 = typename std::iterator_traits<Iter2>::value_type;
	using V1 = std::conditional_t<Vectorize, xsimd::simd_type<T1>, T1>;
	using V2 = std::conditional_t<Vectorize, xsimd::simd_type<T2>, T2>;

	
	[[maybe_unused]] auto carry = make_compensation_carry<OutV, multiplies_result_t<V1, V2>>(reduceOp, init);

	ptrdiff_t idx = 0;
	if (count & 1) {
		const auto v1 = V1(*first1) * uniform_load_unaligned<V2>(std::addressof(*first2));
		std::advance(first1, 1);
		std::advance(first2, -1);

		init += v1;
		idx += 1;
	}
	if (count & 2) {
		const auto v1 = V1(*first1) * uniform_load_unaligned<V2>(std::addressof(*first2));
		std::advance(first1, 1);
		std::advance(first2, -1);
		const auto v2 = V1(*first1) * uniform_load_unaligned<V2>(std::addressof(*first2));
		std::advance(first1, 1);
		std::advance(first2, -1);

		init += v1 + v2;
		idx += 2;
	}
	for (; idx < count; idx += 4) {
		const auto v1 = V1(*first1) * uniform_load_unaligned<V2>(std::addressof(*first2));
		std::advance(first1, 1);
		std::advance(first2, -1);
		const auto v2 = V1(*first1) * uniform_load_unaligned<V2>(std::addressof(*first2));
		std::advance(first1, 1);
		std::advance(first2, -1);
		const auto v3 = V1(*first1) * uniform_load_unaligned<V2>(std::addressof(*first2));
		std::advance(first1, 1);
		std::advance(first2, -1);
		const auto v4 = V1(*first1) * uniform_load_unaligned<V2>(std::addressof(*first2));
		std::advance(first1, 1);
		std::advance(first2, -1);

		if constexpr (is_operator_compensated_v<ReduceOp>) {
			init = reduceOp(carry, init, (v1 + v2) + (v3 + v4));
		}
		else {
			init = reduceOp(init, (v1 + v2) + (v3 + v4));
		}
	}

	return init;
}

template <class Iter1, class Iter2, class IterOut, class ReduceOp = plus_compensated<>>
void ConvolutionReduceVec(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, IterOut firstOut, IterOut lastOut, ptrdiff_t n, bool accumulate = false, ReduceOp reduceOp = plus_compensated<>{}) {
	using T1 = typename std::iterator_traits<Iter1>::value_type;
	using T2 = typename std::iterator_traits<Iter2>::value_type;
	using OutT = typename std::iterator_traits<IterOut>::value_type;

	constexpr bool isVectorized = is_convolution_reduce_vectorized<T1, T2, OutT>::value;
	constexpr ptrdiff_t vectorWidth = isVectorized ? xsimd::simd_traits<OutT>::size : 1;
	using OutV = std::conditional_t<isVectorized, xsimd::simd_type<OutT>, OutT>;

	const ptrdiff_t len1 = std::distance(first1, last1);
	const ptrdiff_t len2 = std::distance(first2, last2);

	// It's better to have input #2 to be longer because then there will be less padding overall.
	if (len2 < len1) {
		return ConvolutionReduceVec(first2, last2, first1, last1, firstOut, lastOut, n, accumulate, reduceOp);
	}

	std::array<T2, vectorWidth * 4 - 2> padding;
	std::fill(padding.begin(), padding.end(), T2(0));
	std::copy(first2, first2 + std::min(vectorWidth, len2), padding.begin() + vectorWidth - 1);
	std::copy(std::reverse_iterator{ last2 }, std::reverse_iterator{ last2 } + std::min(vectorWidth, len2), padding.rbegin() + vectorWidth - 1);

	while (firstOut < lastOut) {
		const ptrdiff_t iterationWidth = std::min(ptrdiff_t(lastOut - firstOut), vectorWidth);
		OutV accumulator = accumulate ? uniform_load_partial_front<OutV>(std::addressof(*firstOut), iterationWidth) : OutV(OutT(0));

		const ptrdiff_t mFirst = std::max(ptrdiff_t(0), n - len2 + 1);
		const ptrdiff_t mLast = std::min(len1, n + vectorWidth);

		const ptrdiff_t mLastPre = std::max(mFirst, std::min(mLast, n + vectorWidth - len2));
		const ptrdiff_t mLastMid = std::max(mLastPre, std::min(n, mLast));
		const ptrdiff_t mLastPost = std::max(mLastMid, mLast);

		const ptrdiff_t paddingPreOffset = ptrdiff_t(padding.size()) - vectorWidth + 1 + n - mFirst - len2;
		const ptrdiff_t midOffset = std::max(ptrdiff_t(0), n - mLastPre);
		const ptrdiff_t paddingPostOffset = n - mLastMid + vectorWidth - 1;
		accumulator = ConvolutionReduceLoop<isVectorized>(first1 + mFirst, padding.data() + paddingPreOffset, accumulator, mLastPre - mFirst, reduceOp);
		accumulator = ConvolutionReduceLoop<isVectorized>(first1 + mLastPre, first2 + midOffset, accumulator, mLastMid - mLastPre, reduceOp);
		accumulator = ConvolutionReduceLoop<isVectorized>(first1 + mLastMid, padding.data() + paddingPostOffset, accumulator, mLastPost - mLastMid, reduceOp);

		uniform_store_partial_front(std::addressof(*firstOut), accumulator, iterationWidth);
		n += iterationWidth;
		firstOut += iterationWidth;
	}
}

} // namespace dspbb::kernels