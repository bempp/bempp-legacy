# Boost
set(BOOST_LIBRARY_DIR "" CACHE STRING "Directory containing Boost Libraries")
set(BOOST_INCLUDE_DIR "" CACHE STRING "Directory containing Boost include files")

# BLAS
set( BLAS_LIBRARIES "" CACHE STRING "Space separated list of full paths to BLAS libs")
#set(ARMADILLO_CMAKE "" CACHE STRING "Directory of Armadillo Cmake file")
find_package( Armadillo REQUIRED REQUIRED)

# Dune

find_library(LIB_DUNE_COMMON common ${CMAKE_SOURCE_DIR}/contrib/dune/dune-common/dune/common/.libs)
find_library(LIB_DUNE_GRID grid ${CMAKE_SOURCE_DIR}/contrib/dune/dune-grid/dune/grid/.libs)

include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-common)
include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-grid)
include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-localfunctions)
include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-foamgrid)
include_directories(${CMAKE_SOURCE_DIR}/contrib/ahmed/Include)
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${BOOST_INCLUDE_DIR})
