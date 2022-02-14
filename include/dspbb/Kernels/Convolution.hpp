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

	template <class U1, class U2, std::enable_if_t<!is_simd_type_v<U1> || !is_simd_type_v<U2>, int> = 0>
	constexpr auto operator()(U1&& accumulator, U2&& increase) const
		-> decltype(math_functions::fma(std::declval<U2>(), std::declval<T>(), std::declval<U1>())) {
		return math_functions::fma(increase, multiplier, accumulator);
	}

	template <class U1, class U2, std::enable_if_t<is_simd_type_v<U1> && is_simd_type_v<U2>, int> = 0>
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

	if (std::min(len1, len2) > 1024) {
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

template <bool Vectorize, class Iter1, class Iter2, class OutV>
inline OutV ConvolutionReduceLoop(Iter1 first1, Iter2 first2, OutV init, ptrdiff_t count) {
	using T1 = typename std::iterator_traits<Iter1>::value_type;
	using T2 = typename std::iterator_traits<Iter2>::value_type;
	using V1 = std::conditional_t<Vectorize, xsimd::simd_type<T1>, T1>;
	using V2 = std::conditional_t<Vectorize, xsimd::simd_type<T2>, T2>;

	for (ptrdiff_t idx = 0; idx != count; ++idx, ++first1, --first2) {
		init += V1{ *first1 } * Load(reinterpret_cast<const V2*>(std::addressof(*first2)));
	}

	return init;
}

template <class Iter1, class Iter2, class IterOut>
void ConvolutionReduceVec(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, IterOut firstOut, IterOut lastOut, ptrdiff_t n, bool accumulate = false) {
	assert(accumulate == false); // Accumulate is not implemented in this function yet, so better make sure it's not misused.

	using T1 = typename std::iterator_traits<Iter1>::value_type;
	using T2 = typename std::iterator_traits<Iter2>::value_type;
	using OutT = typename std::iterator_traits<IterOut>::value_type;

	constexpr bool isVectorized = is_convolution_reduce_vectorized<T1, T2, OutT>::value;
	constexpr ptrdiff_t vectorWidth = isVectorized ? xsimd::simd_traits<OutT>::size : 1;
	using V1 = std::conditional_t<isVectorized, xsimd::simd_type<T1>, T1>;
	using V2 = std::conditional_t<isVectorized, xsimd::simd_type<T2>, T2>;
	using OutV = std::conditional_t<isVectorized, xsimd::simd_type<OutT>, OutT>;

	const ptrdiff_t len1 = std::distance(first1, last1);
	const ptrdiff_t len2 = std::distance(first2, last2);
	std::array<T2, vectorWidth * 4 - 2> padding;
	std::fill(padding.begin(), padding.end(), T2(0));
	std::copy(first2, first2 + std::min(vectorWidth, len2), padding.begin() + vectorWidth - 1);
	std::copy(std::reverse_iterator{ last2 }, std::reverse_iterator{ last2 } + std::min(vectorWidth, len2), padding.rbegin() + vectorWidth - 1);

	for (; firstOut != lastOut; n += vectorWidth) {
		OutV accumulator{ OutT(0) };

		std::fill(accumulator.begin(), accumulator.end(), OutT(0));

		const ptrdiff_t mFirst = std::max(ptrdiff_t(0), n - len2 + 1);
		const ptrdiff_t mLast = std::min(len1, n + vectorWidth);

		ptrdiff_t m = mFirst;
		const ptrdiff_t mLastPre = std::min(mLast, n + vectorWidth - len2);
		const ptrdiff_t mLastMid = std::min(n, mLast);
		const ptrdiff_t mLastPost = mLast;

		for (; m < mLastPre; ++m) {
			const auto access = intptr_t(padding.size()) - vectorWidth + 1 + n - m - len2;
			accumulator = accumulator + V1{ first1[m] } * Load(reinterpret_cast<const V2*>(padding.data() + access));
		}
		for (; m < mLastMid; ++m) {
			const auto access = n - m;
			accumulator = accumulator + V1{ first1[m] } * Load(reinterpret_cast<const V2*>(std::addressof(*(first2 + access))));
		}
		for (; m < mLastPost; ++m) {
			const auto access = n - m + vectorWidth - 1;
			accumulator = accumulator + V1{ first1[m] } * Load(reinterpret_cast<const V2*>(padding.data() + access));
		}

		if (std::distance(firstOut, lastOut) >= vectorWidth) {
			xsimd::store_unaligned(std::addressof(*firstOut), accumulator);
			std::advance(firstOut, vectorWidth);
		}
		else {
			alignas(alignof(OutV)) std::array<OutT, vectorWidth> staging;
			xsimd::store_aligned(staging.data(), accumulator);
			auto stagingIt = staging.begin();
			while (firstOut != lastOut) {
				*firstOut++ = *stagingIt++;
			}
		}
	}
}

} // namespace dspbb::kernels