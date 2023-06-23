add_executable(Benchmark)

target_sources(Benchmark PRIVATE "main.cpp")
target_sources(Benchmark
    PRIVATE
		"Bench_Convolution.cpp"
        "Bench_VectorizedAlgorithms.cpp"
        "Bench_ApplyFilter.cpp"
)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR}/benchmark)
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/benchmark)

find_package(Celero REQUIRED)
target_link_libraries(Benchmark celero DSPBB)