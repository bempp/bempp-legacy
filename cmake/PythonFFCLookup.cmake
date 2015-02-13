# Checks if a python package exists, otherwise installs it

# Adds packages to a python directory
include(PythonPackageLookup)
include(CMakeParseArguments)
include(PythonPackage)
include(PackageLookup)


function(lookup_ffc_package package)
    cmake_parse_arguments(lpp "QUIET;REQUIRED" "VERSION;PATH" "" ${ARGN})
    _python_executable(LOCALPYTHON ${lpp_UNPARSED_ARGUMENTS})
    if(NOT lpp_PATH)
        set(lpp_PATH "${EXTERNAL_ROOT}/python")
    endif()
    set(arguments "")
    set(package_install_name "${package}")
    if(lpp_VERSION)
        list(APPEND arguments VERSION ${lpp_VERSION})
        set(package_install_name "${package}==${lpp_VERSION}")
    endif()
    if(lpp_QUIET)
        list(APPEND arguments QUIET)
    endif()
    find_python_package(${package} ${arguments})
    if(${package}_FOUND)
        return()
    elseif(NOT lpp_QUIET)
        message(STATUS "Will now attempt to install ${package} locally")
    endif()

    # check we have setuptools
    find_python_package(setuptools QUIET ${lpp_UNPARSED_ARGUMENTS})
    if(NOT setuptools_FOUND)
        if(lpp_REQUIRED)
            message(FATAL_ERROR "setuptools not available, cannot install ${package}")
        elseif(NOT lpp_QUIET)
            message(STATUS "setuptools not available, cannot install ${package}")
        endif()
        return()
    endif()

    get_filename_component(lpp_PATH "${lpp_PATH}" ABSOLUTE)
    get_filename_component(EXTERNAL_ROOT "${EXTERNAL_ROOT}" ABSOLUTE)
    if(NOT EXISTS "${lpp_PATH}")
        file(MAKE_DIRECTORY "${lpp_PATH}")
    endif()
    _lpp_check_is_syspath("${LOCALPYTHON}" "${lpp_PATH}")
    file(WRITE "${EXTERNAL_ROOT}/install_${package}.py"
        "from os import environ\n"
        "if 'PYTHONPATH' not in environ:\n"
        "    environ['PYTHONPATH'] = '${lpp_PATH}'\n"
        "elif len(environ['PYTHONPATH']) == 0:\n"
        "    environ['PYTHONPATH'] = '${lpp_PATH}'\n"
        "else:\n"
        "    environ['PYTHONPATH'] += ':${lpp_PATH}'\n"
        "from sys import path, exit\n"
        "path.append('${lpp_PATH}')\n"
        "from setuptools.command.easy_install import main as install\n"
        "result = install(['--install-dir', '${lpp_PATH}', '${package_install_name}'])\n"
        "exit(0 if result == True or result is None else 1)\n"
    )

if(${CMAKE_SYSTEM_NAME} MATCHES "DARWIN")
    file(WRITE "${EXTERNAL_ROOT}/install_ffc.sh"
        "/bin/bash -c \"CXXFLAGS='-stdlib=libc++' MACOSX_DEPLOYMENT_TARGET=10.9 ${LOCALPYTHON} ${EXTERNAL_ROOT}/install_${package}.py\"")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    file(WRITE "${EXTERNAL_ROOT}/install_ffc.sh"
        "/bin/bash -c \"${LOCALPYTHON} ${EXTERNAL_ROOT}/install_${package}.py\"")
endif()
    execute_process(
        COMMAND sh ${EXTERNAL_ROOT}/install_ffc.sh
        RESULT_VARIABLE result
        ERROR_VARIABLE error
        OUTPUT_VARIABLE output
    )
    if(result EQUAL 0)
        if(lpp_REQUIRED)
            list(APPEND arguments REQUIRED)
        endif()
        find_python_package(${package} ${arguments} ${lpp_UNPARSED_ARGUMENTS})
    elseif(lpp_REQUIRED)
        message("output: ${output}\n")
        message("error: ${error}\n")
        message(FATAL_ERROR "Could not install ${package}")
    else()
        message(WARNING "Could not install ${package}")
    endif()
endfunction()
