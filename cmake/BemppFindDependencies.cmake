# Boost
set(BOOST_LIBRARY_DIR "" CACHE PATH "Directory containing Boost Libraries")
set(BOOST_INCLUDE_DIR "" CACHE PATH "Directory containing Boost include files")
include_directories(${BOOST_INCLUDE_DIR})

# BLAS
set(BLAS_LIBRARIES "" CACHE STRING "Space-separated list of full paths to BLAS libs")
#set(ARMADILLO_CMAKE "" CACHE STRING "Directory of Armadillo Cmake file")
find_package( Armadillo REQUIRED REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Ahmed (optional, used only if WITH_AHMED is set)
set(AHMED_INCLUDE_DIR "" CACHE PATH "Full path to the AHMED include directory")
set(AHMED_LIBRARY "" CACHE FILEPATH "Full path to the AHMED library")
set(METIS_LIBRARY "" CACHE FILEPATH "Full path to the METIS library")

# Dune
find_library(LIB_DUNE_COMMON common ${CMAKE_SOURCE_DIR}/contrib/dune/dune-common/dune/common/.libs)
find_library(LIB_DUNE_GRID grid ${CMAKE_SOURCE_DIR}/contrib/dune/dune-grid/dune/grid/.libs)
include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-common)
include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-grid)
include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-localfunctions)
include_directories(${CMAKE_SOURCE_DIR}/contrib/dune/dune-foamgrid)
