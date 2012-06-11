#!/usr/bin/sh
. ../../options.cfg

# mkdir $Main_prefix/bempp/contrib/boost
mkdir -p $Main_prefix/bempp/contrib
mkdir build

cd build

cmake \
 -D CMAKE_CXX_COMPILER:STRING:STRING=$Main_cxx \
 -D CMAKE_C_COMPILER:STRING=$Main_cc \
 -D CMAKE_INSTALL_PREFIX:PATH=$Main_prefix/bempp/contrib/boost \
 -D BUILD_PROJECTS:STRING=test  \
 -D CMAKE_BUILD_TYPE:STRING=Release \
 ..
make
make install
