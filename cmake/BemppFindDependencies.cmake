include(FindPkgConfig)  # use pkg-config to find stuff
include(PythonPackage)  # check for existence of python packages
include(PackageLookup)  # check for existence, or install external projects
include(EnvironmentScript) # scripts to run stuff from build directory

# Creates script for running python with the bempp package available
# Also makes python packages and selected directories available to the build system
add_to_python_path("${PROJECT_BINARY_DIR}/python")
add_to_python_path("${EXTERNAL_ROOT}/python")
add_python_eggs("${PROJECT_SOURCE_DIR}"
    EXCLUDE
        "${PROJECT_SOURCE_DIR}/bempp*egg"
        "${PROJECT_SOURCE_DIR}/Bempp*egg"
)
set(LOCAL_PYTHON_EXECUTABLE "${PROJECT_BINARY_DIR}/localpython.sh")
create_environment_script(
    EXECUTABLE "${PYTHON_EXECUTABLE}"
    PATH "${LOCAL_PYTHON_EXECUTABLE}"
    PYTHON
)

# First, find general packages
find_package(Doxygen)
# Look for lapack and blas
include(lapack_and_blas)
# Look for python libraries corresponding to the python interpreter
# This step is likely not compatible with (automatic) cross-compilation
find_package(CoherentPython REQUIRED)
find_package(Sphinx)
if (WITH_CUDA)
   find_package(CUDA)
endif ()

list(INSERT CMAKE_LOOKUP_PATH 0 ${PROJECT_SOURCE_DIR}/cmake/lookups)
if(WITH_ALUGRID)
    lookup_package(ALUGrid)
endif()
lookup_package(CAIRO REQUIRED)
lookup_package(Boost 1.55 COMPONENTS unit_test_framework REQUIRED)
lookup_package(Armadillo REQUIRED ARGUMENTS TIMEOUT 60)
# ARMA_DONT_USE_WRAPPER means we don't need to include armadillo library
add_definitions(-DARMA_DONT_USE_WRAPPER)
lookup_package(TBB REQUIRED)
lookup_package(Dune REQUIRED COMPONENTS geometry grid localfunctions foamgrid)
lookup_package(SWIG 2.0.4 REQUIRED)
if (SWIG_FOUND AND SWIG_VERSION VERSION_LESS 2.0.7)
    message(WARNING "Swig version 2.0.7 or higher is strongly "
        "recommended to compile BEM++ Python wrappers; "
        "older versions may produce incorrect docstrings"
    )
endif()

# Using cmake_policy does not seem to work here.
set(CMAKE_POLICY_DEFAULT_CMP0012 NEW CACHE STRING "Avoids anoying messages")
unset(arguments)
if(PYPACKED)
    set(arguments ARGUMENTS PYPACKED)
endif()
# Trilinos depends on SWIG, Boost and TBB, so those packages must be looked up
# first.
lookup_package(Trilinos
    DOWNLOAD_BY_DEFAULT REQUIRED CHECK_EXTERNAL
    ${arguments}
)
unset(CMAKE_POLICY_DEFAULT_CMP0012 CACHE)

if("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    # Patchelf allows us to patch the rpath at install time.
    # That's the only time when we know for sure where libraries will be going.
    lookup_package(Patchelf REQUIRED)
endif()

# Then look for python related packages

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
if(NOT EXISTS "${PROJECT_BINARY_DIR}/include/FC.h")
    file(WRITE "${PROJECT_BINARY_DIR}/include/FC.h" "// fake Fortran-C file")
endif()

# Now include all dependency directories once and for all
set(BEMPP_INCLUDE_DIRS
   "${PROJECT_BINARY_DIR}/include/"
   "${PROJECT_BINARY_DIR}/include/bempp"
   ${dune-common_INCLUDE_DIRS}
   ${Trilinos_INCLUDE_DIRS}
   ${Trilinos_TPL_INCLUDE_DIRS}
   ${CAIRO_INCLUDE_DIRS}
)
foreach(component Boost BLAS LAPACK ARMADILLO TBB ALUGrid)
    if(${component}_INCLUDE_DIR)
        list(APPEND BEMPP_INCLUDE_DIRS ${${component}_INCLUDE_DIR})
    endif()
endforeach()


# Add locations of different libraries
add_to_ld_path(
    "${EXTERNAL_ROOT}/lib"
    ${BLAS_LIBRARIES}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE}
    ${TBB_LIBRARY}
    ${TBB_LIBRARY_DEBUG}
    ${TBB_MALLOC_LIBRARY}
    ${TBB_MALLOC_LIBRARY_DEBUG}
    ${CAIRO_LIBRARIES}
)

if(WITH_TESTS)
    include(PythonPackageLookup)
    lookup_python_package(pytest REQUIRED PATH "${EXTERNAL_ROOT}/python")
endif()

# Now adds commands to install external packages
if(EXISTS "${EXTERNAL_ROOT}/lib")
    install(DIRECTORY "${EXTERNAL_ROOT}/lib/"
        DESTINATION "${LIBRARY_INSTALL_PATH}")
endif()
if(EXISTS "${EXTERNAL_ROOT}/include")
    install(DIRECTORY "${EXTERNAL_ROOT}/include/"
        DESTINATION "${INCLUDE_INSTALL_PATH}")
endif()
if(EXISTS "${EXTERNAL_ROOT}/share")
    install(DIRECTORY "${EXTERNAL_ROOT}/share/"
        DESTINATION "${SHARE_INSTALL_PATH}")
endif()
# Trilinos is installed in its own subdirectory, since it is a python package.
if(EXISTS "${EXTERNAL_ROOT}/python/PyTrilinos")
    install(DIRECTORY "${EXTERNAL_ROOT}/python/PyTrilinos"
        DESTINATION "${PYTHON_PKG_DIR}")
endif()

if(WITH_OPENCL)
    find_package(OPENCL REQUIRED)
    list(APPEND BEMPP_INCLUDE_DIRS ${OPENCL_INCLUDE_DIR})
endif()

list(REMOVE_DUPLICATES BEMPP_INCLUDE_DIRS)
include_directories(${BEMPP_INCLUDE_DIRS})
