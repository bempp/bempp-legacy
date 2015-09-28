include(FindPkgConfig)  # use pkg-config to find stuff
include(PackageLookup)  # check for existence, or install external projects
include(EnvironmentScript) # scripts to run stuff from build directory
include(PythonPackageLookup) # adds python packages if not found

# Creates script for running python with the bempp package available
# Also makes python packages and selected directories available to the build
# system
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
find_package(ZLIB REQUIRED)
find_package(Doxygen)
# Look for python libraries corresponding to the python interpreter
# This step is likely not compatible with (automatic) cross-compilation
find_package(CoherentPython REQUIRED)
# Get Python home directory
include(CallPython)
call_python(PYTHON_HOME "import sys; print('%s:%s' % (sys.base_prefix, sys.base_exec_prefix)")

find_package(Sphinx)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(BOOST_MIN_VER 1.57)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(BOOST_MIN_VER 1.54) 
else()
  message(FATAL_ERROR "Windows installation not supported.")
endif()

list(INSERT CMAKE_LOOKUP_PATH 0 ${PROJECT_SOURCE_DIR}/cmake/lookups)
# lookup_package(CAIRO REQUIRED)
lookup_package(Eigen3 REQUIRED)
lookup_package(Boost ${BOOST_MIN_VER} COMPONENTS unit_test_framework filesystem
               program_options system thread iostreams REQUIRED)
lookup_package(TBB REQUIRED)
lookup_package(Dune REQUIRED COMPONENTS geometry grid localfunctions devel )
if (WITH_ALUGRID)
    lookup_package(dune-alugrid REQUIRED)
else()
    lookup_package(dune-foamgrid REQUIRED)
endif()
include("${PROJECT_SOURCE_DIR}/cmake/Dune/local.cmake")

# Using cmake_policy does not seem to work here.
set(CMAKE_POLICY_DEFAULT_CMP0012 NEW CACHE STRING "Avoids anoying messages")
unset(arguments)
if(PYPACKED)
    set(arguments ARGUMENTS PYPACKED)
endif()
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

# Mako is used to generate some of the python bindings
lookup_python_package(mako REQUIRED)
find_program(mako_SCRIPT mako-render HINTS "${EXTERNAL_ROOT}/python")
# Logic for mako should go into this directory
add_to_python_path("${PROJECT_SOURCE_DIR}/python/templates")


# Adds fake FC.h file cos dune incorrectly includes it in dune_config.h
if(NOT EXISTS "${PROJECT_BINARY_DIR}/include/FC.h")
    file(WRITE "${PROJECT_BINARY_DIR}/include/FC.h" "// fake Fortran-C file")
endif()


# Now include all dependency directories once and for all
set(BEMPP_INCLUDE_DIRS
   "${PROJECT_BINARY_DIR}/include/"
   "${PROJECT_BINARY_DIR}/include/bempp"
   ${dune-common_INCLUDE_DIRS}
   ${PYTHON_INCLUDE_DIR}
   ${NUMPY_INCLUDE_DIRS}
   ${dune-alugrid_INCLUDE_DIRS}
   ${dune-foamgrid_INCLUDE_DIRS}
)

foreach(component Boost TBB EIGEN3)
    if(${component}_INCLUDE_DIR)
        list(APPEND BEMPP_INCLUDE_DIRS ${${component}_INCLUDE_DIR})
    endif()
endforeach()


# Add locations of different libraries
add_to_ld_path(
    ${ZLIB_LIBRARIES}
    "${EXTERNAL_ROOT}/lib"
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE}
    ${TBB_LIBRARY}
    ${TBB_LIBRARY_DEBUG}
    ${TBB_MALLOC_LIBRARY}
    ${TBB_MALLOC_LIBRARY_DEBUG}
)

lookup_python_package(Cython VERSION 0.21 REQUIRED PATH "${EXTERNAL_ROOT}/python")
if(WITH_TESTS)
    include(AddPyTest)
    setup_pytest("${EXTERNAL_ROOT}/python" "${PROJECT_BINARY_DIR}/py.test.sh")
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
        DESTINATION "${SHARE_INSTALL_PATH}"
        PATTERN "doc" EXCLUDE)
endif()


list(REMOVE_DUPLICATES BEMPP_INCLUDE_DIRS)
include_directories(${BEMPP_INCLUDE_DIRS})
