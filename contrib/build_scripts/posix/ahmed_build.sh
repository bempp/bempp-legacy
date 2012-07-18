#!/usr/bin/sh
set -e
. ../../options.cfg

mkdir -p build

cd build
cmake \
 -D CMAKE_CXX_COMPILER:STRING:STRING=$Main_cxx \
 -D CMAKE_C_COMPILER:STRING=$Main_cc \
 -D CMAKE_INSTALL_PREFIX:PATH=$Main_prefix/bempp/contrib/ahmed \
 -D CMAKE_BUILD_TYPE:STRING=Release \
 -D CMAKE_CXX_FLAGS:STRING=$Main_cxxflags \
 -D CMAKE_C_FLAGS:STRING=$Main_cflags \
 -D BLAS_LIBS:STRING=$BLAS_lib \
 -D LAPACK_LIBS:STRING=$LAPACK_lib \
 -D ENABLE_METIS:BOOL=OFF \
 ..
make
make install
