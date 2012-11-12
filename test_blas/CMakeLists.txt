cmake_minimum_required(VERSION 2.8)
project(test_blas)

# BLAS
set(BLAS_LIBRARIES "" CACHE STRING "Semicolon-separated list of full paths to BLAS libs")

add_executable(test_blas test_blas.cpp)
add_executable(test_sdot_returns_float test_sdot_returns_float.cpp)
add_executable(test_sdot_returns_double test_sdot_returns_double.cpp)
add_executable(test_zdotc test_zdotc.cpp)
add_executable(test_zdotc_g77 test_zdotc_g77.cpp)

target_link_libraries(test_blas ${BLAS_LIBRARIES})
target_link_libraries(test_sdot_returns_float ${BLAS_LIBRARIES})
target_link_libraries(test_sdot_returns_double ${BLAS_LIBRARIES})
target_link_libraries(test_zdotc ${BLAS_LIBRARIES})
target_link_libraries(test_zdotc_g77 ${BLAS_LIBRARIES})

