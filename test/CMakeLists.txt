add_executable(UnitTest)

target_sources(UnitTest PRIVATE main.cpp)

if (${CMAKE_CXX_COMPILER} MATCHES "MSVC")
	target_compile_options(UnitTest PRIVATE "/bigobj")
endif()

target_sources(UnitTest 
	PRIVATE
		"Filtering/FIR/Test_Descs.cpp"
		"Filtering/IIR/Test_BandTransforms.cpp"
		"Filtering/IIR/Test_Descs.cpp"
		"Filtering/IIR/Test_Realizations.cpp"
		"Filtering/Test_FIR.cpp"
		"Filtering/Test_IIR.cpp"
		"Filtering/Test_MeasureFilter.cpp"
		"Filtering/Test_Polyphase.cpp"
		"Filtering/Test_Resample.cpp"
		"Filtering/Test_Windowing.cpp"
		"Generators/Test_Generators.cpp"
		"Kernels/Test_Convolution.cpp" 
		"Kernels/Test_Numeric.cpp" 
		"Kernels/Test_Numeric.cpp"
		"LTISystems/Test_DiscretizationTransforms.cpp"
		"LTISystems/Test_Systems.cpp"
		"Math/Test_Convolution.cpp"
		"Math/Test_EllipticFunctions.cpp"
		"Math/Test_FFT.cpp"
		"Math/Test_Functions.cpp"
		"Math/Test_OverlapAdd.cpp"
		"Math/Test_Polynomials.cpp"
		"Math/Test_Rational.cpp"
		"Math/Test_RootTransforms.cpp"
		"Math/Test_Solvers.cpp"
		"Math/Test_Statistics.cpp"
		"Primitives/Test_Signal.cpp"
		"Primitives/Test_SignalArithmetic.cpp"
		"Primitives/Test_SignalView.cpp"
		"Utility/Test_Interval.cpp"
)

find_package(Catch2 REQUIRED)

target_link_libraries(UnitTest DSPBB)
target_link_libraries(UnitTest Catch2::Catch2)