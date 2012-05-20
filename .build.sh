#!/usr/bin/sh
source ./options.cfg

mkdir build
cd build
cmake \
    -D CMAKE_INSTALL_PREFIX:STRING=$Main_prefix \
    -D CMAKE_BUILD_TYPE:STRING=Release \
    -D BOOST_INCLUDE_DIR:STRING=$Boost_include_dir \
    -D BOOST_UNIT_TEST_LIB:STRING=$Boost_unit_test_lib \
    -D BLAS_LIBRARIES:STRING=$BLAS_lib \
    -D LAPACK_LIBRARIES:STRING=$LAPACK_lib \
    -D ARMADILLO_INCLUDE_DIR:STRING=$Armadillo_include_dir \
    -D TBB_INCLUDE_DIR:STRING=$Tbb_include_dir \
    -D TBB_LIBRARY:STRING=$Tbb_lib \
    -D TBB_LIBRARY_DEBUG:STRING=$Tbb_lib_debug \
    -D WITH_AHMED:BOOL=ON \
    -D AHMED_INCLUDE_DIR:STRING=$AHMED_include_dir \
    -D AHMED_LIB:STRING=$AHMED_lib \
    -D METIS_LIB:STRING=$AHMED_metis_lib \
    -D TRILINOS_CMAKE_PATH:STRING=$Trilinos_cmake_path \
..
