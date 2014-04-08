include(FindPkgConfig)

# First, find general packages
find_package(Doxygen)
find_package(Sphinx)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(CoherentPython REQUIRED)
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


# Then look for python related packages
include(PythonPackage)
function(find_or_fail package what)
    find_python_package(${package})
    if(NOT ${package}_FOUND)
        message("*********")
        message("${package} is required to ${what}")
        message("It can likely be installed with pip")
        message("*********")
        message(FATAL_ERROR "Aborting")
    endif()
endfunction()

# first looks for python package, second for linkage/include stuff
find_or_fail(numpy "by Purify's python bindings")
find_package(Numpy REQUIRED)

lookup_package(SWIG 2.0.4 REQUIRED)
if (SWIG_VERSION VERSION_LESS 2.0.7)
    message(WARNING "Swig version 2.0.7 or higher is strongly "
        "recommended to compile BEM++ Python wrappers; "
        "older versions may produce incorrect docstrings"
    )
endif ()


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
    ${Boost_INCLUDE_DIR}
    ${TBB_INCLUDE_DIR}
    ${dune-common_INCLUDE_DIRS}
    ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS}
)
if(ARMADILLO_INCLUDE_DIR)
    include_directories(${ARMADILLO_INCLUDE_DIR})
endif()


if(WITH_TESTS)
    # Adds a virtual environment
    find_or_fail(virtualenv "to run the unit-tests for the python bindings")
    include(PythonVirtualEnv)
    # Add paths so we can run tests effectively
    add_to_python_path("${PROJECT_BINARY_DIR}/python")
    add_to_ld_path("${EXTERNAL_ROOT}/lib")
    add_package_to_virtualenv(pytest)
endif()
