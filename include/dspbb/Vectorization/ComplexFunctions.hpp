#pragma once

#pragma warning(push)
#pragma warning(disable: 4800)
#include <xsimd/xsimd.hpp>
#pragma warning(pop)

namespace dspbb {
namespace complex_functions {

	template <class T>
	struct AbsVec {
		static constexpr size_t stride = xsimd::simd_traits<std::complex<T>>::size;
		using complex_vector = xsimd::batch<std::complex<T>, stride>;
		using real_vector = xsimd::batch<T, stride>;
		void operator()(T* out, const std::complex<T>* in) {
			complex_vector vin;
			vin.load_unaligned(in);
			const real_vector vout = xsimd::abs(vin);
			vout.store_unaligned(out);
		}
	};

	template <class T>
	struct Abs {
		void operator()(T* out, const std::complex<T>* in) {
			*out = std::abs(*in);
		}
	};

	
	template <class T>
	struct ArgVec {
		static constexpr size_t stride = xsimd::simd_traits<std::complex<T>>::size;
		using complex_vector = xsimd::batch<std::complex<T>, stride>;
		using real_vector = xsimd::batch<T, stride>;
		void operator()(T* out, const std::complex<T>* in) {
			complex_vector vin;
			vin.load_unaligned(in);
			const real_vector vout = xsimd::arg(vin);
			vout.store_unaligned(out);
		}
	};

	template <class T>
	struct Arg {
		void operator()(T* out, const std::complex<T>* in) {
			*out = std::arg(*in);
		}
	};

		
	template <class T>
	struct RealVec {
		static constexpr size_t stride = xsimd::simd_traits<std::complex<T>>::size;
		using complex_vector = xsimd::batch<std::complex<T>, stride>;
		using real_vector = xsimd::batch<T, stride>;
		void operator()(T* out, const std::complex<T>* in) {
			complex_vector vin;
			vin.load_unaligned(in);
			const real_vector vout = xsimd::real(vin);
			vout.store_unaligned(out);
		}
	};

	template <class T>
	struct Real {
		void operator()(T* out, const std::complex<T>* in) {
			*out = std::real(*in);
		}
	};
	
	template <class T>
	struct ImagVec {
		static constexpr size_t stride = xsimd::simd_traits<std::complex<T>>::size;
		using complex_vector = xsimd::batch<std::complex<T>, stride>;
		using real_vector = xsimd::batch<T, stride>;
		void operator()(T* out, const std::complex<T>* in) {
			complex_vector vin;
			vin.load_unaligned(in);
			const real_vector vout = xsimd::imag(vin);
			vout.store_unaligned(out);
		}
	};

	template <class T>
	struct Imag {
		void operator()(T* out, const std::complex<T>* in) {
			*out = std::imag(*in);
		}
	};

} // namespace complex_functions
} // namespace dspbb