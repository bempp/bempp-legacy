#!/bin/sh
. ./.options.cfg

export LD_LIBRARY_PATH=$Main_prefix/bempp/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$Main_prefix/bempp/lib:$DYLD_LIBRARY_PATH

mkdir $Bempp_build_dir
cd $Bempp_build_dir
CXXFLAGS="$Main_cxxflags" CFLAGS="$Main_cflags" cmake \
    -D CMAKE_CXX_COMPILER:STRING=$Main_cxx \
    -D CMAKE_C_COMPILER:STRING=$Main_cc \
    -D CMAKE_INSTALL_PREFIX:STRING=$Main_prefix \
    -D CMAKE_BUILD_TYPE:STRING=$Bempp_build_type \
    -D WITH_AHMED:STRING=$AHMED_with_ahmed \
    -D BOOST_INCLUDE_DIR:STRING=$Boost_include_dir \
    -D BOOST_UNIT_TEST_LIB:STRING=$Boost_unit_test_lib \
    -D BLAS_LIBRARIES:STRING=$BLAS_lib \
    -D LAPACK_LIBRARIES:STRING=$LAPACK_lib \
    -D ARMADILLO_INCLUDE_DIR:STRING=$Armadillo_include_dir \
    -D TBB_INCLUDE_DIR:STRING=$Tbb_include_dir \
    -D TBB_LIBRARY:STRING=$Tbb_lib \
    -D TBB_LIBRARY_DEBUG:STRING=$Tbb_lib_debug \
    -D AHMED_INCLUDE_DIR:STRING=$AHMED_include_dir \
    -D AHMED_LIB:STRING=$AHMED_lib \
    -D PYTHON_EXECUTABLE:STRING=$Python_exe \
    -D PYTHON_INCLUDE_DIR:STRING=$Python_include_dir \
    -D PYTHON_NUMPY_INCLUDE_DIR:STRING=$Python_numpy_include_dir \
    -D PYTHON_LIBRARY:STRING=$Python_lib \
    -D TRILINOS_CMAKE_PATH:STRING=$Trilinos_cmake_path \
$Main_root_dir
