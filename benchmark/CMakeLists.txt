add_executable(Benchmark)

target_sources(Benchmark PRIVATE "main.cpp")
target_sources(Benchmark
    PRIVATE
		"Bench_Convolution.cpp"
        "Bench_VectorizedAlgorithms.cpp"
)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR}/benchmark)
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/benchmark)

conan_cmake_configure(
	REQUIRES
		celero/2.8.2
	GENERATORS cmake_find_package
)
conan_cmake_install(PATH_OR_REFERENCE . BUILD missing SETTINGS ${conan_settings} OUTPUT_QUIET)


find_package(Celero REQUIRED)
target_link_libraries(Benchmark celero::celero DSPBB)