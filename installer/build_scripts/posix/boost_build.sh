#!/bin/sh
set -e
. ../../.options.cfg

mkdir -p build

cd build

CXXFLAGS="$Boost_cxxflags" CFLAGS="$Boost_cflags" "$CMake_exe" \
 -D CMAKE_CXX_COMPILER:STRING:STRING="$Boost_cxx" \
 -D CMAKE_C_COMPILER:STRING="$Boost_cc" \
 -D CMAKE_INSTALL_PREFIX:PATH="$Main_prefix/bempp" \
 -D BUILD_PROJECTS:STRING=test  \
 -D CMAKE_BUILD_TYPE:STRING=Release \
 -D BOOST_INCLUDE_INSTALL_DIR:STRING=include \
 -D BOOST_LIB_INSTALL_DIR:STRING=lib \
 ..
