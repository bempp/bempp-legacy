list(INSERT CMAKE_LOOKUP_PATH 0 ${PROJECT_SOURCE_DIR}/cmake/lookups)
lookup_package(Boost COMPONENTS unit_test_framework REQUIRED)
lookup_package(Armadillo REQUIRED)
lookup_package(TBB REQUIRED)
lookup_package(Dune REQUIRED DOWNLOAD_BY_DEFAULT
    COMPONENTS geometry grid localfunctions foamgrid)
lookup_package(Trilinos
    REQUIRED
    ARGUMENTS
        URL /Users/mdavezac/workspace/bempp/trilinos-11.6.1-Source.tar.bz2
        MD5 b97d882535fd1856599b1c7338f5b45a
)
lookup_package(SWIG 2.0.4 REQUIRED)

if(Boost_FOUND)
  if(CMAKE_BUILD_TYPE STREQUAL "Release" OR MAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(BOOST_UNIT_TEST_LIB ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE})
  else()
    set(BOOST_UNIT_TEST_LIB ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG})
  endif()

  set(
    BOOST_UNIT_TEST_LIB ${BOOST_UNIT_TEST_LIB}
    CACHE INTERNAL
    "Path to unit test framework"
  )
  set(
    BOOST_INCLUDE_DIR ${Boost_INCLUDE_DIR}
    CACHE INTERNAL
    "Path to boost include directory"
  )
endif()

# BLAS
find_package(BLAS REQUIRED)
# find_package blas does not look for cblas.h
if(NOT BLAS_INCLUDE_DIR)
    find_path(BLAS_INCLUDE_DIR cblas.h)
endif()

# LAPACK
find_package(LAPACK REQUIRED)
if(NOT LAPACK_INCLUDE_DIR)
    find_path(LAPACK_INCLUDE_DIR lapack.h)
endif()

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
file(GLOB_RECURSE DUNE_HEADERS ${CMAKE_INSTALL_PREFIX}/bempp/include/dune/*.hh)
