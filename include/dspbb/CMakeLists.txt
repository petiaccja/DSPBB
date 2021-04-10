# Enumerate sources
add_library(DSPBB INTERFACE)

message("${CMAKE_CURRENT_SOURCE_DIR}/..")
target_include_directories(DSPBB INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/..")

find_package(xsimd REQUIRED)
target_link_libraries(DSPBB INTERFACE xsimd::xsimd)