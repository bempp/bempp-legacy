# Checks for python package
# See https://github.com/UCL/GreatCMakeCookOff/wiki for information

# First check for python executable
include(FindPackageHandleStandardArgs)
include(Utilities)
include(CMakeParseArguments)

find_package(PythonInterp REQUIRED)

function(_python_executable OUTVAR)
    cmake_parse_arguments(PYEXEC "LOCAL" "PYTHON_EXECUTABLE" "" ${ARGN})
    if(PYEXEC_LOCAL AND PYEXEC_PYTHON_EXECUTABLE)
        message(FATAL_ERROR "Cannot use LOCAL and PYTHON arguments together")
    endif()
    if(PYEXEC_LOCAL AND NOT LOCAL_PYTHON_EXECUTABLE)
        message(FATAL_ERROR "PythonVirtualEnv not included yet.")
    endif()
    if(PYEXEC_LOCAL)
        set(${OUTVAR} ${LOCAL_PYTHON_EXECUTABLE} PARENT_SCOPE)
    elseif(PYEXEC_PYTHON_EXECUTABLE)
        set(${OUTVAR} ${PYEXEC_PYTHON_EXECUTABLE} PARENT_SCOPE)
    elseif(NOT LOCAL_PYTHON_EXECUTABLE)
        set(${OUTVAR} ${PYTHON_EXECUTABLE} PARENT_SCOPE)
    elseif(EXISTS "${LOCAL_PYTHON_EXECUTABLE}")
        set(${OUTVAR} ${LOCAL_PYTHON_EXECUTABLE} PARENT_SCOPE)
    endif()
    set(${OUTVAR}_UNPARSED_ARGUMENTS ${PYEXEC_UNPARSED_ARGUMENTS} PARENT_SCOPE)
endfunction()

function(find_python_package PACKAGE)
    string(TOUPPER "${PACKAGE}" PACKAGE_UPPER)
    if(${PACKAGE_UPPER}_FOUND)
        set(${PACKAGE}_FOUND ${${PACKAGE_UPPER}_FOUND} PARENT_SCOPE)
        return()
    endif()
    cmake_parse_arguments(${PACKAGE}_FIND
        "REQUIRED;EXACT" "VERSION" "" ${ARGN})
    cmake_parse_arguments(PYPACK
        "QUIET" "WORKING_DIRECTORY;ERROR_MESSAGE" ""
        ${${PACKAGE}_FIND_UNPARSED_ARGUMENTS})
    _python_executable(LOCALPYTHON ${PYPACK_UNPARSED_ARGUMENTS})
    if(NOT PYPACK_WORKING_DIRECTORY)
        set(PYPACK_WORKING_DIRECTORY "${PROJECT_BINARY_DIR}")
    endif()
    if(NOT PYPACK_ERROR_MESSAGE)
        set(PYPACK_ERROR_MESSAGE "Python module ${PACKAGE} could not be found.")
    endif()
    if(PYPACK_QUIET)
        set(${PACKAGE}_FIND_QUIETLY TRUE)
        set(${PACKAGE_UPPER}_FIND_QUIETLY TRUE)
    else()
        set(${PACKAGE}_FIND_QUIETLY FALSE)
        set(${PACKAGE_UPPER}_FIND_QUIETLY FALSE)
    endif()
    # Unset prior variables
    unset(${PACKAGE}_LOCATION CACHE)
    unset(${PACKAGE}_FOUND CACHE)
    unset(${PACKAGE_UPPER}_FOUND CACHE)

    execute_process(
        COMMAND ${LOCALPYTHON} -c
            "import ${PACKAGE};print(getattr(${PACKAGE}, '__version__', ''))"
        WORKING_DIRECTORY "${PYPACK_WORKING_DIRECTORY}"
        RESULT_VARIABLE PACKAGE_WAS_FOUND
        ERROR_VARIABLE ERROR
        OUTPUT_VARIABLE OUTPUT
    )
    if("${PACKAGE_WAS_FOUND}" STREQUAL "0")
        string(STRIP "${OUTPUT}" string_version)
        set(arguments
            "import ${PACKAGE}"
            "from os.path import dirname"
            "print(dirname(${PACKAGE}.__file__))"
        )
        execute_process(
            COMMAND ${LOCALPYTHON} -c "${arguments}"
            WORKING_DIRECTORY "${PYPACK_WORKING_DIRECTORY}"
            RESULT_VARIABLE LOCATION_WAS_FOUND
            ERROR_VARIABLE ERROR
            OUTPUT_VARIABLE OUTPUT
        )
        if("${LOCATION_WAS_FOUND}" STREQUAL "0")
            string(STRIP "${OUTPUT}" OUTPUT)
            FILE(TO_CMAKE_PATH ${OUTPUT} ${PACKAGE}_LOCATION)
        endif()
    endif()
    find_package_handle_standard_args(${PACKAGE}
        REQUIRED_VARS ${PACKAGE}_LOCATION
        VERSION_VAR string_version
        FAIL_MESSAGE "${PYPACK_ERROR_MESSAGE}"
    )
    if(NOT "${PACKAGE}" STREQUAL "${PACKAGE_UPPER}")
        set(${PACKAGE}_FOUND ${${PACKAGE_UPPER}_FOUND} PARENT_SCOPE)
    endif()
    set(${PACKAGE_UPPER}_FOUND "${${PACKAGE_UPPER}_FOUND}" CACHE INTERNAL "")
    if(${PACKAGE_UPPER}_FOUND)
        set(${PACKAGE}_LOCATION "${${PACKAGE}_LOCATION}"
            CACHE PATH "Location of ${PACKAGE}")
        set(${PACKAGE}_VERSION_STRING "${string_version}"
            CACHE STRING "Version of ${PACKAGE}")
    else()
        set(${PACKAGE}_LOCATION "${${PACKAGE}_LOCATION}"
            CACHE PATH "Location of ${PACKAGE}")
    endif()
endfunction()
