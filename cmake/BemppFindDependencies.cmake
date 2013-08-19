# Boost
set(BOOST_INCLUDE_DIR "" CACHE PATH "Directory containing Boost header files")
set(BOOST_UNIT_TEST_LIB "" CACHE PATH "Full path to Boost unit test library")

# BLAS
set(BLAS_LIBRARIES "" CACHE STRING "Semicolon-separated list of full paths to BLAS libs")
set(BLAS_INCLUDE_DIR "" CACHE PATH "Directory containing BLAS header files (used only with MKL, GotoBLAS and OpenBLAS)")

# LAPACK
set(LAPACK_LIBRARIES "" CACHE STRING "Semicolon-separated list of full paths to LAPACK libs")
set(LAPACK_INCLUDE_DIR "" CACHE PATH "Directory containing LAPACK header files (used only with MKL, GotoBLAS and OpenBLAS)")

# ARMADILLO
set(ARMADILLO_INCLUDE_DIR "" CACHE PATH "Directory containing Armadillo header files")

# Threading building blocks
set(TBB_INCLUDE_DIR "" CACHE PATH "Directory containing Intel TBB header files")
set(TBB_LIBRARY "" CACHE PATH "Full path to the TBB library")
set(TBB_LIBRARY_DEBUG "" CACHE PATH "Full path to the TBB debug library")
set(TBB_MALLOC_LIBRARY "" CACHE PATH "Full path to the TBB malloc library")
set(TBB_MALLOC_LIBRARY_DEBUG "" CACHE PATH "Full path to the TBB malloc debug library")

# Ahmed (optional, used only if WITH_AHMED is set)
if (WITH_AHMED)
    set(AHMED_INCLUDE_DIR "" CACHE PATH "Full path to the AHMED include directory")
    set(AHMED_LIB "" CACHE PATH "Full path to AHMED library")
endif ()

# CUDA support
if (WITH_CUDA)
   FIND_PACKAGE(CUDA)
endif ()

# Dune
find_library(LIB_DUNE_COMMON dunecommon ${CMAKE_INSTALL_PREFIX}/bempp/lib)
find_library(LIB_DUNE_GRID dunegrid ${CMAKE_INSTALL_PREFIX}/bempp/lib)
file(GLOB_RECURSE DUNE_HEADERS ${CMAKE_INSTALL_PREFIX}/bempp/include/dune/*.hh)

# Trilinos
set(TRILINOS_CMAKE_PATH "" CACHE PATH "Directory containing TrilinosConfig.cmake")
include(${TRILINOS_CMAKE_PATH}/TrilinosConfig.cmake)


