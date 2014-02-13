# simple CMake file to install tbb from its extracted download file to places we want.
cmake_minimum_required(VERSION 2.8)

# Find patch files and apply
file(GLOB HEADER_FILES dune/foamgrid/*.hh)
install(FILES ${HEADER_FILES} DESTINATION include/dune/foamgrid)

file(GLOB HEADER_FILES dune/foamgrid/foamgrid/*.hh)
install(FILES ${HEADER_FILES} DESTINATION include/dune/foamgrid/foamgrid)
