# Enumerate sources
add_library(DSPBB INTERFACE)

message("${CMAKE_CURRENT_SOURCE_DIR}/..")
target_include_directories(DSPBB INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/..")

find_package(xsimd REQUIRED)
find_package(Eigen3 REQUIRED)
target_link_libraries(DSPBB INTERFACE xsimd::xsimd)
target_link_libraries(DSPBB INTERFACE Eigen3::Eigen3)