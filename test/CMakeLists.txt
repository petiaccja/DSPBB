# Enumerate sources
add_executable(UnitTest)

target_sources(UnitTest PRIVATE main.cpp)

target_sources(UnitTest 
	PRIVATE
		"ComputeKernels/Test_Kernels.cpp" 
		"Filtering/FIR/Test_FIRDescs.cpp"
		"Filtering/IIR/Test_BandTransforms.cpp"
		"Filtering/IIR/Test_Descs.cpp"
		"Filtering/IIR/Test_Kernels.cpp"
		"Filtering/IIR/Test_Realizations.cpp"
		"Filtering/Test_FilterParameters.cpp"
		"Filtering/Test_FirFilter.cpp"
		"Filtering/Test_IIR.cpp"
		"Filtering/Test_Interpolation.cpp"
		"Filtering/Test_OverlapAdd.cpp"
		"Filtering/Test_Polyphase.cpp"
		"Filtering/Test_Windowing.cpp"
		"Generators/Test_Generators.cpp"
		"LTISystems/Test_DiscretizationTransforms.cpp"
		"Math/Test_Convolution.cpp"
		"Math/Test_FFT.cpp"
		"Math/Test_Functions.cpp"
		"Math/Test_Polynomials.cpp"
		"Math/Test_RootTransforms.cpp"
		"Math/Test_Statistics.cpp"
		"Primitives/Test_Signal.cpp"
		"Primitives/Test_SignalArithmetic.cpp"
		"Primitives/Test_SignalView.cpp"
)

find_package(Catch2 REQUIRED)

target_link_libraries(UnitTest DSPBB)
target_link_libraries(UnitTest Catch2::Catch2)