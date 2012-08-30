cmake_minimum_required(VERSION 2.8)
project(test_blas)

# BLAS
set(BLAS_LIBRARIES "" CACHE STRING "Semicolon-separated list of full paths to BLAS libs")
set(LAPACK_LIBRARIES "" CACHE STRING "Semicolon-separated list of full paths to LAPACK libs")

add_executable(test_lapack test_lapack.cpp)

target_link_libraries(test_lapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

