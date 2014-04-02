list(INSERT CMAKE_LOOKUP_PATH 0 ${PROJECT_SOURCE_DIR}/cmake/lookups)
lookup_package(Boost COMPONENTS unit_test_framework REQUIRED)
lookup_package(Armadillo REQUIRED DOWNLOAD_WARNING)
# lookup_package(TBB REQUIRED)
# lookup_package(Dune)
# lookup_package(Trillinos
#     ARGUMENTS LOCATION /Users/mdavezac/workspace/bempp
# )
lookup_package(SWIG 2.0.4 REQUIRED DOWNLOAD_WARNING)
return()

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
set(BLAS_LIBRARIES "" CACHE STRING "Semicolon-separated list of full paths to BLAS libs")
set(BLAS_INCLUDE_DIR "" CACHE PATH "Directory containing BLAS header files (used only with MKL, GotoBLAS and OpenBLAS)")

# LAPACK
set(LAPACK_LIBRARIES "" CACHE STRING "Semicolon-separated list of full paths to LAPACK libs")
set(LAPACK_INCLUDE_DIR "" CACHE PATH "Directory containing LAPACK header files (used only with MKL, GotoBLAS and OpenBLAS)")

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

