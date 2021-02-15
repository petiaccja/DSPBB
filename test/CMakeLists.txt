# Enumerate sources
add_executable(UnitTest ${sources} ${headers})

target_sources(UnitTest PRIVATE main.cpp "Test_Functions.cpp")

target_sources(UnitTest 
	PRIVATE
		Test_Convolution.cpp
		Test_FFT.cpp
		Test_FirFilter.cpp
		Test_Interpolation.cpp
		Test_Polyphase.cpp
		Test_Signal.cpp
		Test_SignalView.cpp
		Test_WindowFunctions.cpp
		Test_MathFunctions.cpp
		Test_Algorithm.cpp
 "Test_Statistics.cpp")


target_include_directories(UnitTest PUBLIC ${CMAKE_SOURCE_DIR}/support)
target_include_directories(UnitTest PUBLIC ${CMAKE_SOURCE_DIR}/include)