#pragma once

#ifdef _MSC_VER
	#pragma warning(push)
	#pragma warning(disable : 4800 4244)
#endif
#include <xsimd/xsimd.hpp>
#ifdef _MSC_VER
	#pragma warning(pop)
#endif

namespace dspbb::kernels {

template <class T, class U>
inline T uniform_load_unaligned(const U* mem) {
	if constexpr (xsimd::is_batch<std::decay_t<T>>::value) {
		return T::load_unaligned(mem);
	}
	else {
		return *mem;
	}
}

template <class T, class U>
inline void uniform_store_unaligned(U* mem, const T& value) {
	if constexpr (xsimd::is_batch<std::decay_t<T>>::value) {
		value.store_unaligned(mem);
	}
	else {
		*mem = value;
	}
}

template <class VecT, class T>
inline VecT uniform_load_partial_front(const T* data, size_t count) {
	if constexpr (!xsimd::is_batch<std::decay_t<VecT>>::value) {
		return *data;
	}
	else {
		constexpr auto vectorWidth = xsimd::revert_simd_traits<std::decay_t<VecT>>::size;
		if (count == vectorWidth) {
			return VecT::load_unaligned(data);
		}
		std::array<T, vectorWidth> extended;
		std::copy(data, data + count, extended.begin());
		return VecT::load_unaligned(extended.data());
	}
}

template <class VecT, class T>
inline void uniform_store_partial_front(T* data, const VecT& v, size_t count) {
	if constexpr (!xsimd::is_batch<std::decay_t<VecT>>::value) {
		*data = v;
	}
	else {
		constexpr auto vectorWidth = xsimd::revert_simd_traits<std::decay_t<VecT>>::size;
		if (count == vectorWidth) {
			v.store_unaligned(data);
			return;
		}
		alignas(alignof(VecT)) std::array<T, vectorWidth> extended;
		v.store_unaligned(extended.data());
		std::copy(extended.begin(), extended.begin() + count, data);
	}
}


} // namespace dspbb::kernels