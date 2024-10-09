add_library(DSPBB INTERFACE)

target_sources(DSPBB
	INTERFACE FILE_SET headers TYPE HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/.." FILES
        # Filtering
        "FilterUtility.hpp"
        "FIR.hpp"
        "IIR.hpp"
        "MeasureFilter.hpp"
        "Polyphase.hpp"
        "Resample.hpp"
        "Windowing.hpp"
        
        # FIR
        "FIR/BandTransforms.hpp"
        "FIR/Descs.hpp"
        "FIR/Filter.hpp"
        "FIR/LeastSquares.hpp"
        "FIR/Windowed.hpp"
        
        # IIR
        "IIR/BandTransforms.hpp"
        "IIR/Butterworth.hpp"
        "IIR/Chebyshev.hpp"
        "IIR/Descs.hpp"
        "IIR/Elliptic.hpp"
        "IIR/Filter.hpp"
        "IIR/Realizations.hpp"
        
        # Generators
        "Spaces.hpp"
        "Waveforms.hpp"
        
        # Kernels
        "Convolution.hpp"
        "Functors.hpp"
        "Math.hpp"
        "Numeric.hpp"
        "Utility.hpp"
        
        # LTISystems
        "DiscretizationTransforms.hpp"
        "Systems.hpp"
        
        # Math
        "Convolution.hpp"
        "DotProduct.hpp"
        "EllipticFunctions.hpp"
        "FFT.hpp"
        "Functions.hpp"
        "OverlapAdd.hpp"
        "Polynomials.hpp"
        "Rational.hpp"
        "RootTransforms.hpp"
        "Solvers.hpp"
        "Statistics.hpp"
        
        # Primitives
        "Signal.hpp"
        "SignalArithmetic.hpp"
        "Traits.hpp"
        "SignalView.hpp"
        
        # Utility
        "Interval.hpp"
        "Numbers.hpp"
        "TypeTraits.hpp"
)

target_compile_features(DSPBB INTERFACE cxx_std_20)
set_target_properties(DSPBB PROPERTIES VERIFY_INTERFACE_HEADER_SETS ON INTERFACE_HEADER_SETS_TO_VERIFY headers)

find_package(xsimd REQUIRED)
find_package(pocketfft REQUIRED)
find_package(Eigen3 REQUIRED)
target_link_libraries(DSPBB INTERFACE xsimd)
target_link_libraries(DSPBB INTERFACE pocketfft::pocketfft)
target_link_libraries(DSPBB INTERFACE Eigen3::Eigen)