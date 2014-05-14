# simple CMake file to install tbb from its extracted download file to places we want.
cmake_minimum_required(VERSION 2.8)

# Find patch files and apply
file(GLOB HEADER_FILES dune/foamgrid/*.hh)
install(FILES ${HEADER_FILES} DESTINATION include/dune/foamgrid)

file(GLOB HEADER_FILES dune/foamgrid/foamgrid/*.hh)
install(FILES ${HEADER_FILES} DESTINATION include/dune/foamgrid/foamgrid)

file(WRITE "dune-foamgrid.pc"
   "prefix=${CMAKE_INSTALL_PREFIX}\n"
   "exec_prefix=\${prefix}\n"
   "libdir=\${exec_prefix}/lib\n"
   "includedir=\${prefix}/include\n"
   "CXX=${CMAKE_CXX_FLAGS}\n"
   "CC=${CMAKE_C_FLAGS}\n\n"
   "DEPENDENCIES= dune-common   dune-grid  \n"
   "Name: dune-foamgrid\n"
   "Version: 2.1\n"
   "Description: dune-foamgrid module\n"
   "URL: http://dune-project.org/\n"
   "Requires: dune-common dune-grid\n"
   "Libs: \n"
   "Cflags: \n"
)
install(FILES dune-foamgrid.pc DESTINATION lib/pkgconfig)
