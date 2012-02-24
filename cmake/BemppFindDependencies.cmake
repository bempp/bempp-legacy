# Boost
find_package( Boost REQUIRED)

# BLAS
set( BLAS_LIBRARIES "" CACHE STRING "Space separated list of full paths to BLAS libs")

# Dune

find_library(LIB_DUNE_COMMON common ${CMAKE_SOURCE_DIR}/contrib/dune/dune-common/dune/common/.libs)
find_library(LIB_DUNE_GRID grid ${CMAKE_SOURCE_DIR}/contrib/dune/dune-grid/dune/grid/.libs)

include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-common)
include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-grid)
include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-localfunctions)
include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-foamgrid)
include_directories(${CMAKE_SOURCE_DIR}/contrib/ahmed/Include)

