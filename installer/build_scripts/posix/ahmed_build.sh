#!/bin/sh
set -e
. ../../.options.cfg

mkdir -p build

export DYLD_LIBRARY_PATH="$Main_prefix/bempp/lib:$DYLD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$Main_prefix/bempp/lib:$LD_LIBRARY_PATH"

cd build
CXXFLAGS="$AHMED_cxxflags" CFLAGS="$AHMED_cflags" "$CMake_exe" \
 -D CMAKE_CXX_COMPILER:STRING:STRING="$AHMED_cxx" \
 -D CMAKE_C_COMPILER:STRING="$AHMED_cc" \
 -D CMAKE_INSTALL_PREFIX:PATH="$Main_prefix/bempp" \
 -D CMAKE_BUILD_TYPE:STRING=Release \
 -D BLAS_LIBS:STRING=$BLAS_lib \
 -D LAPACK_LIBS:STRING=$LAPACK_lib \
 -D ENABLE_METIS:BOOL=OFF \
 -D ENABLE_64BIT:BOOL=$AHMED_enable64 \
 ..
