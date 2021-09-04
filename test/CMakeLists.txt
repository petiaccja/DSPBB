# Enumerate sources
add_executable(UnitTest)

target_sources(UnitTest PRIVATE main.cpp)

target_sources(UnitTest 
	PRIVATE
		"Filtering/Test_Convolution.cpp"
		"Filtering/Test_FFT.cpp"
		"Filtering/Test_FirFilter.cpp"
		"Generators/Test_Generators.cpp"
		"Filtering/Test_Interpolation.cpp"
		"Filtering/Test_OverlapAdd.cpp"
		"Filtering/Test_Polyphase.cpp"
		"Primitives/Test_Signal.cpp"
		"Primitives/Test_SignalView.cpp"
		"Primitives/Test_SignalArithmetic.cpp"
		"Filtering/Test_WindowFunctions.cpp"		
		
		
		"Math/Test_Statistics.cpp"
		"Math/Test_Functions.cpp"
		"Vectorization/Test_Kernels.cpp" 
		)

find_package(Catch2 REQUIRED)

target_link_libraries(UnitTest DSPBB)
target_link_libraries(UnitTest Catch2::Catch2)