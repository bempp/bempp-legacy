#!/bin/sh
set -e
. ../../options.cfg

mkdir -p build

cd build

CXXFLAGS="$Main_cxxflags" CFLAGS="$Main_cflags" cmake \
 -D CMAKE_CXX_COMPILER:STRING:STRING=$Main_cxx \
 -D CMAKE_C_COMPILER:STRING=$Main_cc \
 -D CMAKE_INSTALL_PREFIX:PATH=$Main_prefix/bempp \
 -D BUILD_PROJECTS:STRING=test  \
 -D CMAKE_BUILD_TYPE:STRING=Release \
 -D BOOST_INCLUDE_INSTALL_DIR:STRING=include \
 -D BOOST_LIB_INSTALL_DIR:STRING=lib \
 ..
