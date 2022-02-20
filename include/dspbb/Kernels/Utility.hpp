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

template <class T>
struct is_simd_type : std::false_type {};

template <class Scalar, auto... Params>
struct is_simd_type<xsimd::batch<Scalar, Params...>> : std::true_type {};

template <class T>
constexpr bool is_simd_type_v = is_simd_type<T>::value;

} // namespace dspbb::kernels