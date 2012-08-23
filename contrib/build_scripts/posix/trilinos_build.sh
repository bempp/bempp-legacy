#!/usr/bin/sh
set -e
. ../../.options.cfg

# mkdir $Main_prefix/bempp/contrib/trilinos
mkdir -p $Main_prefix/bempp/contrib
mkdir -p build
cd build
CXXFLAGS="$Main_cxxflags" CFLAGS="$Main_cflags" cmake \
    -D CMAKE_CXX_COMPILER:STRING=$Main_cxx \
    -D CMAKE_C_COMPILER:STRING=$Main_cc \
    -D CMAKE_BUILD_TYPE:STRING=Release \
    -D CMAKE_INSTALL_PREFIX:PATH=$Main_prefix/bempp \
    -D BUILD_SHARED_LIBS:BOOL=ON \
    -D TPL_TBB_LIBRARIES:STRING=$Tbb_lib \
    -D TPL_TBB_INCLUDE_DIRS:PATH=$Tbb_include_dir \
    -D TPL_LAPACK_LIBRARIES:STRING=$LAPACK_lib \
    -D TPL_BLAS_LIBRARIES:STRING=$BLAS_lib \
    -D TPL_ENABLE_BLAS:BOOL=ON \
    -D TPL_ENABLE_LAPACK:BOOL=ON \
    -D TPL_ENABLE_TBB:BOOL=ON \
    -D TPL_ENABLE_Boost:BOOL=ON \
    -D TPL_Boost_INCLUDE_DIRS:STRING=$Boost_include_dir \
    -D Trilinos_ENABLE_Amesos:BOOL=ON \
    -D Trilinos_ENABLE_Belos:BOOL=ON \
    -D Trilinos_ENABLE_Epetra:BOOL=ON \
    -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
    -D Trilinos_ENABLE_RTOp:BOOL=ON \
    -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
    -D Trilinos_ENABLE_Teuchos:BOOL=ON \
    -D Trilinos_ENABLE_Thyra:BOOL=ON \
    -D Trilinos_ENABLE_ThyraCore:BOOL=ON \
    -D Trilinos_ENABLE_ThyraEpetraAdapters:BOOL=ON \
    -D Trilinos_ENABLE_ThyraEpetraExtAdapters:BOOL=ON \
    -D Trilinos_ENABLE_ThyraTpetraAdapters:BOOL=ON \
    -D Trilinos_ENABLE_Tpetra:BOOL=ON \
    -D Trilinos_ENABLE_Fortran:BOOL=OFF \
    -D RTOp_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
    -D Stratimikos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
    -D Thyra_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
    -D Tpetra_INST_COMPLEX_FLOAT:BOOL=ON \
    -D Tpetra_INST_FLOAT:BOOL=ON \
..

