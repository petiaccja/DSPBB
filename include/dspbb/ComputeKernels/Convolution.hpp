#pragma once

#include "VectorizedAlgorithms.hpp"

#include <utility>


namespace dspbb::kernels {


template <class R, class U, class V>
void Convolution(R* out, const U* u, const V* v, size_t lenU, size_t lenV, size_t first, size_t count) {
	if (lenU < lenV) {
		return Convolution(out, v, u, lenV, lenU, first, count);
	}

	memset(out, 0, sizeof(R) * count);
	for (size_t i = 0; i < lenV; ++i) {
		const auto scale = v[i];
		const intptr_t uoffset = std::max(intptr_t(0), intptr_t(first) - intptr_t(i));
		const intptr_t ooffset = std::max(intptr_t(0), intptr_t(i) - intptr_t(first));
		const intptr_t ccount = std::max(intptr_t(0), std::min(intptr_t(count) - ooffset, intptr_t(lenU) - uoffset));
		BinaryOperationVectorized(out + ooffset, out + ooffset, u + uoffset, ccount, [scale](const auto& a, const auto& b) { return a + b * scale; });
	}
}



} // namespace dspbb::kernels