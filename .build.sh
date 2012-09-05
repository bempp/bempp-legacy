#!/bin/sh
. ./.options.cfg

export LD_LIBRARY_PATH="$Main_prefix/bempp/lib:$LD_LIBRARY_PATH"
export DYLD_LIBRARY_PATH="$Main_prefix/bempp/lib:$DYLD_LIBRARY_PATH"

mkdir "$Bempp_build_dir"
cd $Bempp_build_dir
CXXFLAGS="$Bempp_cxxflags" CFLAGS="$Bempp_cflags" cmake \
    -D CMAKE_CXX_COMPILER:STRING="$Bempp_cxx" \
    -D CMAKE_C_COMPILER:STRING="$Bempp_cc" \
    -D CMAKE_INSTALL_PREFIX:PATH="$Main_prefix" \
    -D CMAKE_BUILD_TYPE:STRING="$Bempp_build_type" \
    -D WITH_AHMED:STRING="$AHMED_with_ahmed" \
    -D BOOST_INCLUDE_DIR:PATH="$Boost_include_dir" \
    -D BOOST_UNIT_TEST_LIB:PATH="$Boost_unit_test_lib" \
    -D BLAS_LIBRARIES:STRING="$BLAS_lib" \
    -D LAPACK_LIBRARIES:STRING="$LAPACK_lib" \
    -D ARMADILLO_INCLUDE_DIR:PATH="$Armadillo_include_dir" \
    -D TBB_INCLUDE_DIR:PATH="$Tbb_include_dir" \
    -D TBB_LIBRARY:PATH="$Tbb_lib" \
    -D TBB_LIBRARY_DEBUG:PATH="$Tbb_lib_debug" \
    -D AHMED_INCLUDE_DIR:PATH="$AHMED_include_dir" \
    -D AHMED_LIB:PATH="$AHMED_lib" \
    -D PYTHON_EXECUTABLE:PATH="$Python_exe" \
    -D PYTHON_INCLUDE_DIR:PATH="$Python_include_dir" \
    -D PYTHON_NUMPY_INCLUDE_DIR:PATH="$Python_numpy_include_dir" \
    -D PYTHON_LIBRARY:PATH="$Python_lib" \
    -D SWIG_EXECUTABLE:PATH="$Swig_exe" \
    -D TRILINOS_CMAKE_PATH:PATH="$Trilinos_cmake_path" \
    -D WITH_MKL:BOOL="$MKL_enable_mkl" \
"$Main_root_dir"
