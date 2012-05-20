# Boost
set(BOOST_INCLUDE_DIR "" CACHE STRING "Directory containing Boost include files")
set(BOOST_UNIT_TEST_LIB "" CACHE STRING "Full path to Boost unit test library")

# BLAS
set(BLAS_LIBRARIES "" CACHE STRING "Semicolon-separated list of full paths to BLAS libs")

# LAPACK
set(LAPACK_LIBRARIES "" CACHE STRING "Semicolon-separated list of full paths to LAPACK libs")


# ARMADILLO

set(ARMADILLO_INCLUDE_DIR "" CACHE STRING "Include directory for Armadillo")



# Threading building blocks

set(TBB_INCLUDE_DIR "" CACHE STRING "Full path to Intel TBB include directory")
set(TBB_LIBRARY "" CACHE STRING "Full path to the TBB library")
set(TBB_LIBRARY_DEBUG "" CACHE STRING "Full path to the TBB debug library")

# Ahmed (optional, used only if WITH_AHMED is set)
if (WITH_AHMED)
	set(AHMED_INCLUDE_DIR "" CACHE STRING "Full path to the AHMED include directory")
	set(AHMED_LIB "" CACHE STRING "Full path to AHMED library")
	set(METIS_LIB "" CACHE STRING "Full path to METIS library (only requiref if AHMED is built with Metis)")
endif ()

# Dune

find_library(LIB_DUNE_COMMON common ${CMAKE_INSTALL_PREFIX}/bempp/contrib/dune/dune-common/dune/common/.libs)
find_library(LIB_DUNE_GRID grid ${CMAKE_INSTALL_PREFIX}/bempp/contrib/dune/dune-grid/dune/grid/.libs)

file(GLOB_RECURSE DUNE_FOAMGRID_HEADERS ${CMAKE_INSTALL_PREFIX}/bempp/contrib/dune/dune-foamgrid/dune/foamgrid/*.hh)
file(GLOB_RECURSE DUNE_GRID_HEADERS ${CMAKE_INSTALL_PREFIX}/bempp/contrib/dune/dune-grid/dune/grid/*.hh)
file(GLOB_RECURSE DUNE_COMMON_HEADERS ${CMAKE_INSTALL_PREFIX}/bempp/contrib/dune/dune-common/dune/common/*.hh)
file(GLOB_RECURSE DUNE_LOCALFUNCTIONS_HEADERS ${CMAKE_INSTALL_PREFIX}/bempp/contrib/dune/dune-localfunctions/dune/localfunctions/*.hh)

# Trilinos

set(TRILINOS_CMAKE_PATH "" CACHE STRING "Directory containing TrilinosConfig.cmake")
include(${TRILINOS_CMAKE_PATH}/TrilinosConfig.cmake)


