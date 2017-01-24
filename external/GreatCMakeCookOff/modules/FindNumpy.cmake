# - Find the NumPy libraries
# This module finds if NumPy is installed, and sets the following variables
# indicating where it is.
#
# TODO: Update to provide the libraries and paths for linking npymath lib.
#
#  NUMPY_FOUND               - was NumPy found
#  NUMPY_VERSION             - the version of NumPy found as a string
#  NUMPY_VERSION_MAJOR       - the major version number of NumPy
#  NUMPY_VERSION_MINOR       - the minor version number of NumPy
#  NUMPY_VERSION_PATCH       - the patch version number of NumPy
#  NUMPY_INCLUDE_DIRS        - path to the NumPy include files

# Modified from script by Continuum Analytics, Inc.
if(NUMPY_INCLUDE_DIRS) # Don't look twice
    string(REGEX REPLACE "\\." ";" _NUMPY_VERSION_LIST ${NUMPY_VERSION})
    list(GET _NUMPY_VERSION_LIST 0 NUMPY_VERSION_MAJOR)
    list(GET _NUMPY_VERSION_LIST 1 NUMPY_VERSION_MINOR)
    list(GET _NUMPY_VERSION_LIST 2 NUMPY_VERSION_PATCH)
    return()
endif()

# Finding NumPy involves calling the Python interpreter
find_package(CoherentPython)
if(NOT PYTHON_EXECUTABLE)
    if(NumPy_FIND_REQUIRED)
        message(FATAL_ERROR "Could not find python interpreter.")
    endif()
    return()
endif()

execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "import numpy as n; print(n.__version__); print(n.get_include());"
    RESULT_VARIABLE _NUMPY_SEARCH_SUCCESS
    OUTPUT_VARIABLE _NUMPY_VALUES_OUTPUT
    ERROR_VARIABLE _NUMPY_ERROR_VALUE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT _NUMPY_SEARCH_SUCCESS MATCHES 0)
    if(NumPy_FIND_REQUIRED)
        message(FATAL_ERROR
            "NumPy import failure:\n${_NUMPY_ERROR_VALUE}")
    endif()
    set(NUMPY_FOUND FALSE)
    return()
endif()

# Convert the process output into a list
string(REGEX REPLACE ";" "\\\\;" _NUMPY_VALUES ${_NUMPY_VALUES_OUTPUT})
string(REGEX REPLACE "\n" ";" _NUMPY_VALUES ${_NUMPY_VALUES})
# Just in case there is unexpected output from the Python command.
list(GET _NUMPY_VALUES -2 NUMPY_VERSION)
list(GET _NUMPY_VALUES -1 NUMPY_INCLUDE_DIRS)

string(REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+" _VER_CHECK "${NUMPY_VERSION}")
if("${_VER_CHECK}" STREQUAL "")
    # The output from Python was unexpected. Raise an error always
    # here, because we found NumPy, but it appears to be corrupted somehow.
    message(FATAL_ERROR
        "Requested version and include path from NumPy, got instead:\n${_NUMPY_VALUES_OUTPUT}\n")
    return()
endif()

# Make sure all directory separators are '/'
string(REGEX REPLACE "\\\\" "/" NUMPY_INCLUDE_DIRS ${NUMPY_INCLUDE_DIRS})

# Get the major and minor version numbers
string(STRIP "${NUMPY_VERSION}" NUMPY_VERSION)
string(REGEX REPLACE "\\." ";" _NUMPY_VERSION_LIST ${NUMPY_VERSION})
list(GET _NUMPY_VERSION_LIST 0 NUMPY_VERSION_MAJOR)
list(GET _NUMPY_VERSION_LIST 1 NUMPY_VERSION_MINOR)
list(GET _NUMPY_VERSION_LIST 2 NUMPY_VERSION_PATCH)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(NumpyLibrary
    REQUIRED_VARS
        NUMPY_INCLUDE_DIRS
        NUMPY_VERSION_MAJOR
        NUMPY_VERSION_MINOR
    VERSION_VAR NUMPY_VERSION
)

if(NUMPYLIBRARY_FOUND)
    set(NUMPYLIBRARY_FOUND TRUE CACHE BOOL "Numpy library was found")
    set(NumpyLibrary_FOUND TRUE CACHE BOOL "Numpy library was found")
    set(NUMPY_INCLUDE_DIRS
        "${NUMPY_INCLUDE_DIRS}" CACHE
        PATH "Path to numpy includes"
    )
    set(NUMPY_VERSION "${NUMPY_VERSION}" CACHE INTERNAL "Numpy version")
else()
    return()
endif()

if(NOT no_numpy_feature_tests)
    ## Now check some features of numpy c api
    macro(numpy_feature_test OUTVARNAME testfilename testname)
        ## try to compile and run
        ## Using Release flags because MSCrapware fails otherwise.
        try_compile(
          ${OUTVARNAME}
          ${CMAKE_BINARY_DIR}
          ${CMAKE_CURRENT_LIST_DIR}/numpy/${testfilename}
          COMPILE_DEFINITIONS -I${PYTHON_INCLUDE_DIRS}  -I${NUMPY_INCLUDE_DIRS}
                              -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
          CMAKE_FLAGS -DLINK_LIBRARIES:STRING=${PYTHON_LIBRARIES}
                      -DCMAKE_CXX_FLAGS_DEBUG:STRING="${CMAKE_CXX_FLAGS_RELEASE}"
                      -DCMAKE_C_FLAGS_DEBUG:STRING="${CMAKE_C_FLAGS_RELEASE}"
                      -DCMAKE_EXE_LINKER_FLAGS_DEBUG:STRING="${CMAKE_EXE_LINKER_FLAGS_RELEASE}"
          OUTPUT_VARIABLE NUMPY_TESTCOMPILE
        )
        ## display results
        if(NOT Numpy_FIND_QUIETLY)
            message (STATUS "[NumPy] ${testname} = ${${OUTVARNAME}}")
        endif()
        set(${OUTVARNAME} ${${OUTVARNAME}}
            CACHE BOOL
            "Numpy feature check: ${testname}"
        )
    endmacro()

    numpy_feature_test(NUMPY_NPY_LONG_DOUBLE test_numpy_long_double.cc "Long double exists")
    numpy_feature_test(NUMPY_NPY_BOOL test_numpy_ubyte.cc "Bool is a separate type")
    numpy_feature_test(NUMPY_NPY_ARRAY test_numpy_is_noarray.c "NPY_ARRAY_* macros exist")
    numpy_feature_test( NUMPY_NPY_ENABLEFLAGS test_numpy_has_enableflags.c
                        "PyArray_ENABLEFLAGS exists" )
endif()
