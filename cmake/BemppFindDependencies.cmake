include(FindPkgConfig)
include(PackageLookup)

# First, find general packages
find_package(Doxygen)
find_package(CBLAS REQUIRED)
include(BlasThreads)
find_package(LAPACK REQUIRED)
find_package(CoherentPython REQUIRED)
find_package(Sphinx)
if (WITH_CUDA)
   find_package(CUDA)
endif ()
if(NOT LAPACK_INCLUDE_DIR)
    find_path(LAPACK_INCLUDE_DIR clapack.h)
endif()

list(INSERT CMAKE_LOOKUP_PATH 0 ${PROJECT_SOURCE_DIR}/cmake/lookups)
if(WITH_ALUGRID)
    lookup_package(ALUGrid)
endif()
lookup_package(Boost 1.35 COMPONENTS unit_test_framework REQUIRED)
lookup_package(Armadillo REQUIRED)
lookup_package(TBB DOWNLOAD_BY_DEFAULT REQUIRED CHECK_EXTERNAL)
lookup_package(Dune REQUIRED DOWNLOAD_BY_DEFAULT
    COMPONENTS geometry grid localfunctions foamgrid
    CHECK_EXTERNAL
)
# Using cmake_policy does not seem to work here.
set(CMAKE_POLICY_DEFAULT_CMP0012 NEW CACHE STRING "Avoids anoying messages")
lookup_package(Trilinos DOWNLOAD_BY_DEFAULT REQUIRED CHECK_EXTERNAL
    ARGUMENTS PYPACKED)
unset(CMAKE_POLICY_DEFAULT_CMP0012 CACHE)

# Then look for python related packages
lookup_package(SWIG 2.0.4 REQUIRED)
if (SWIG_FOUND AND SWIG_VERSION VERSION_LESS 2.0.7)
    message(WARNING "Swig version 2.0.7 or higher is strongly "
        "recommended to compile BEM++ Python wrappers; "
        "older versions may produce incorrect docstrings"
    )
endif()

include(PythonPackage)

# first looks for python package, second for linkage/include stuff
find_python_package(numpy REQUIRED
    ERROR_MESSAGE
        "numpy is required by the BEM++ python bindings"
)
find_package(Numpy REQUIRED)


# Ahmed (optional, used only if WITH_AHMED is set)
if (WITH_AHMED)
    set(AHMED_INCLUDE_DIR "" CACHE PATH "Full path to the AHMED include directory")
    set(AHMED_LIB "" CACHE PATH "Full path to AHMED library")
endif ()

# Adds fake FC.h file cos dune incorrectly includes it in dune_config.h
file(WRITE ${PROJECT_BINARY_DIR}/include/FC.h "// fake Fortran-C file")

# Now include all dependency directories once and for all
include_directories(
    ${PROJECT_BINARY_DIR}/include/
    ${PROJECT_BINARY_DIR}/include/bempp
    ${PROJECT_SOURCE_DIR}/lib
    ${dune-common_INCLUDE_DIRS}
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
)
foreach(component Boost BLAS LAPACK ARMADILLO TBB ALUGrid)
    if(${component}_INCLUDE_DIR)
        include_directories(${${component}_INCLUDE_DIR})
    endif()
endforeach()
add_definitions(-DARMA_USE_LAPACK -DARMA_USE_BLAS)


# Creates script for running python with the bempp package available
include(EnvironmentScript)
add_to_ld_path("${EXTERNAL_ROOT}/lib")
add_to_python_path("${PROJECT_BINARY_DIR}/python")
add_to_python_path("${EXTERNAL_ROOT}/python")
set(LOCAL_PYTHON_EXECUTABLE "${PROJECT_BINARY_DIR}/localpython.sh")
create_environment_script(
    EXECUTABLE "${PYTHON_EXECUTABLE}"
    PATH "${LOCAL_PYTHON_EXECUTABLE}"
    PYTHON
)

if(WITH_TESTS)
    include(PythonPackageLookup)
    lookup_python_package(pytest REQUIRED PATH "${EXTERNAL_ROOT}/python")
endif()

# Now adds commands to install external packages
if(EXISTS "${EXTERNAL_ROOT}/lib")
    install(DIRECTORY "${EXTERNAL_ROOT}/lib/" DESTINATION lib)
endif()
if(EXISTS "${EXTERNAL_ROOT}/include")
    install(DIRECTORY "${EXTERNAL_ROOT}/include/" DESTINATION include)
endif()
if(EXISTS "${EXTERNAL_ROOT}/share")
    install(DIRECTORY "${EXTERNAL_ROOT}/share/" DESTINATION share)
endif()

if(WITH_OPENCL)
    find_package(OPENCL REQUIRED)
    include_directories(${OPENCL_INCLUDE_DIR})
endif()
