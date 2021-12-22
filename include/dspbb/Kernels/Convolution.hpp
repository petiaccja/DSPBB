#pragma once

#include "../Utility/Interval.hpp"
#include "Numeric.hpp"

#include <iterator>
#include <utility>


namespace dspbb::kernels {


template <class R, class U, class V>
void Convolution(R* out, const U* u, const V* v, size_t lenU, size_t lenV, size_t first, size_t count, bool clearOut = true) {
	if (lenU < lenV) {
		return Convolution(out, v, u, lenV, lenU, first, count, clearOut);
	}

	if (clearOut) {
		memset(out, 0, sizeof(R) * count);
	}
	for (size_t i = 0; i < lenV; ++i) {
		const auto scale = v[i];
		const intptr_t uoffset = std::max(intptr_t(0), intptr_t(first) - intptr_t(i));
		const intptr_t ooffset = std::max(intptr_t(0), intptr_t(i) - intptr_t(first));
		const intptr_t ccount = std::max(intptr_t(0), std::min(intptr_t(count) - ooffset, intptr_t(lenU) - uoffset));
		Transform(out + ooffset, out + ooffset + ccount, u + uoffset, out + ooffset,
				  [scale](const auto& a, const auto& b) -> plus_result_t<decltype(a), multiplies_result_t<decltype(b), decltype(scale)>> { return a + b * scale; });
	}
}

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
				  [multiplier](const auto& out, const auto& val2) { return out + val2 * multiplier; });
	}
}

template <class Iter1, class Iter2, class IterOut>
void ConvolutionAccumulate(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, IterOut firstOut, IterOut lastOut, ptrdiff_t n, bool accumulate = false) {
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

template <class Iter1, class Iter2, class IterOut>
void ConvolutionAccumulate_Vec(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, IterOut firstOut, IterOut lastOut, ptrdiff_t n, bool accumulate = false) {
	const ptrdiff_t len1 = std::distance(first1, last1);
	const ptrdiff_t len2 = std::distance(first2, last2);

	constexpr ptrdiff_t vectorWidth = 8;
	using T2 = typename std::iterator_traits<Iter2>::value_type;
	using OutT = typename std::iterator_traits<IterOut>::value_type;

	std::array<T2, vectorWidth * 4 - 2> padding;
	std::fill(padding.begin(), padding.end(), T2(0));
	std::copy(first2, first2 + std::min(vectorWidth, len2), padding.begin() + vectorWidth - 1);
	std::copy(std::reverse_iterator{ last2 }, std::reverse_iterator{ last2 } + std::min(vectorWidth, len2), padding.rbegin() + vectorWidth - 1);

	for (; firstOut < lastOut; n += vectorWidth) {
		std::array<OutT, vectorWidth> accumulator;

		std::fill(accumulator.begin(), accumulator.end(), OutT(0));

		ptrdiff_t mFirst = std::max(ptrdiff_t(0), n - len2 + 1);
		ptrdiff_t mLast = std::min(len1, n + vectorWidth);

		for (; n - mFirst + vectorWidth > len2; ++mFirst) {
			const auto access = intptr_t(padding.size()) - vectorWidth + 1 + n - mFirst - len2;
			for (ptrdiff_t i = 0; i < vectorWidth; ++i) {
				accumulator[i] += first1[mFirst] * padding[access + i];
			}
		}
		for (; mFirst < n && mFirst < mLast; ++mFirst) {
			const auto access = n - mFirst;
			for (ptrdiff_t i = 0; i < vectorWidth; ++i) {
				accumulator[i] += first1[mFirst] * first2[access + i];
			}
		}
		for (; mFirst < mLast; ++mFirst) {
			const auto access = n - mFirst + vectorWidth - 1;
			for (ptrdiff_t i = 0; i < vectorWidth; ++i) {
				accumulator[i] += first1[mFirst] * padding[access + i];
			}
		}

		for (size_t i = 0; i < vectorWidth && firstOut < lastOut; ++i, ++firstOut) {
			*firstOut = accumulator[i];
		}
	}
}

} // namespace dspbb::kernels