#!/usr/bin/sh
source ../../options.cfg

mkdir $Main_prefix/bempp/contrib/boost
mkdir build


cd build

cmake \
 -DCMAKE_INSTALL_PREFIX:PATH=$Main_prefix/bempp/contrib/boost \
 -DBUILD_PROJECTS:STRING=test  \
 -DCMAKE_BUILD_TYPE:STRING=Release \
 ..
make
make install
