include(FindPkgConfig)

# First, find general packages
find_package(Doxygen)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
if (WITH_CUDA)
   find_package(CUDA)
endif ()
if(NOT BLAS_INCLUDE_DIR)
    find_path(BLAS_INCLUDE_DIR cblas.h)
endif()
if(NOT LAPACK_INCLUDE_DIR)
    find_path(LAPACK_INCLUDE_DIR clapack.h)
endif()

list(INSERT CMAKE_LOOKUP_PATH 0 ${PROJECT_SOURCE_DIR}/cmake/lookups)
lookup_package(Boost COMPONENTS unit_test_framework REQUIRED)
lookup_package(Armadillo REQUIRED)
lookup_package(TBB DOWNLOAD_BY_DEFAULT REQUIRED)
lookup_package(Dune REQUIRED DOWNLOAD_BY_DEFAULT
    COMPONENTS geometry grid localfunctions foamgrid
)
lookup_package(Trilinos
    DOWNLOAD_BY_DEFAULT
    REQUIRED
    ARGUMENTS
        URL /Users/mdavezac/workspace/bempp/trilinos-11.6.1-Source.tar.bz2
        MD5 b97d882535fd1856599b1c7338f5b45a
)
lookup_package(SWIG 2.0.4 REQUIRED)


### The following sets bempp specific variables, for ease of use
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


# Ahmed (optional, used only if WITH_AHMED is set)
if (WITH_AHMED)
    set(AHMED_INCLUDE_DIR "" CACHE PATH "Full path to the AHMED include directory")
    set(AHMED_LIB "" CACHE PATH "Full path to AHMED library")
endif ()

# Dune
file(GLOB_RECURSE DUNE_HEADERS ${CMAKE_INSTALL_PREFIX}/bempp/include/dune/*.hh)

# Adds fake FC.h file cos dune incorrectly includes it in dune_config.h
file(WRITE ${PROJECT_BINARY_DIR}/fakery/FC.h "// fake Fortran-C file")

# Now include all dependency directories once and for all
include_directories(
    ${PROJECT_BINARY_DIR}/fakery/
    ${CMAKE_BINARY_DIR}/include
    ${CMAKE_SOURCE_DIR}/lib
    ${BLAS_INCLUDE_DIR}
    ${LAPACK_INCLUDE_DIR}
    ${BOOST_INCLUDE_DIR}
    ${TBB_INCLUDE_DIR}
    ${dune-common_INCLUDE_DIRS}
    ${ARMADILLO_INCLUDE_DIR}
    ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS}
)
