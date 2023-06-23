## DSPBB CMake
cmake_minimum_required(VERSION 3.24.0)

project(DSPBB)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Cpp flags
enable_language(CXX)

option(ENABLE_LLVM_COV "Adds compiler flags to generate LLVM source-based code coverage. Only works with Clang." OFF)
option(ENABLE_LLVM_ADDRESS_SANITIZER "Adds compiler flags to generate LLVM source-based code coverage. Only works with Clang." OFF)
option(ENABLE_LLVM_MEMORY_SANITIZER "Adds compiler flags to generate LLVM source-based code coverage. Only works with Clang." OFF)
option(BUILD_EXAMPLES "Include or exclude examples in the project." ON)
option(BUILD_BENCHMARKS "Include or exclude benchmarks in the project." ON)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
	if (ENABLE_LLVM_COV)
		add_compile_options("-fprofile-instr-generate" "-fcoverage-mapping" "-mllvm" "-enable-name-compression=false")
		add_link_options("-fprofile-instr-generate" "-fcoverage-mapping")
	endif()
	if (ENABLE_LLVM_ADDRESS_SANITIZER)
		add_compile_options("-fsanitize=address")
		add_link_options("-fsanitize=address")
	endif()
	if (ENABLE_LLVM_MEMORY_SANITIZER)
		add_compile_options("-fsanitize=memory")
		add_link_options("-fsanitize=memory")
	endif()
	if ("${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}" STREQUAL "MSVC")
		set(RELWITHDEBINFO_OPTIONS "/Ob2")
		set(RELWITHDEBINFO_REMOVE "/Ob1")
		foreach(option ${RELWITHDEBINFO_REMOVE})
			string(REPLACE ${option} "" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
		endforeach()
		add_compile_options("$<$<CONFIG:RELWITHDEBINFO>:${RELEASE_OPTIONS}>")
	endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
	add_compile_options("/MP")
	set(RELWITHDEBINFO_OPTIONS "/Ob2")
	set(RELWITHDEBINFO_REMOVE "/Ob1")
	foreach(option ${RELWITHDEBINFO_REMOVE})
		string(REPLACE ${option} "" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
	endforeach()
	add_compile_options("$<$<CONFIG:RELWITHDEBINFO>:${RELEASE_OPTIONS}>")
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
link_libraries(Threads::Threads)

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Bin)

# Subdirectories
add_subdirectory(include/dspbb)
add_subdirectory(test)
if (${BUILD_EXAMPLES})
	add_subdirectory(examples)
endif()
if (${BUILD_BENCHMARKS})
	add_subdirectory(benchmark)
endif()

