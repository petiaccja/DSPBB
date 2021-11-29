add_executable(01_Basics "01_Basics.cpp")
add_executable(02_Statistics "02_Statistics.cpp")

if (${CMAKE_CXX_COMPILER} MATCHES "MSVC")
	target_compile_options(01_Basics PRIVATE "/bigobj")
	target_compile_options(02_Statistics PRIVATE "/bigobj")
endif()

target_link_libraries(01_Basics DSPBB)
target_link_libraries(02_Statistics DSPBB)