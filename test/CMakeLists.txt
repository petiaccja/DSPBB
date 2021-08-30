# Enumerate sources
add_executable(UnitTest)

target_sources(UnitTest PRIVATE main.cpp)

target_sources(UnitTest 
	PRIVATE
		Test_Convolution.cpp
		Test_FFT.cpp
		Test_FirFilter.cpp
		Test_Generators.cpp
		Test_Interpolation.cpp
		Test_OverlapAdd.cpp
		Test_Polyphase.cpp
		Test_Signal.cpp
		Test_SignalView.cpp
		Test_SignalArithmetic.cpp
		Test_WindowFunctions.cpp		
		Test_Algorithm.cpp
		Test_Arithmetic.cpp
		Test_Statistics.cpp
		Test_Functions.cpp
		Test_Kernels.cpp 
		Test_Hilbert.cpp)

find_package(Catch2 REQUIRED)

target_link_libraries(UnitTest DSPBB)
target_link_libraries(UnitTest Catch2::Catch2)