# Checks if a python package exists, otherwise installs it

# Adds packages to a python directory
include(PythonPackageLookup)
include(CMakeParseArguments)
include(PythonPackage)
include(PackageLookup)


function(lookup_fiat_package package)
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
        "result = install(['-f','http://www.bempp.org/files/','--install-dir', '${lpp_PATH}', '${package_install_name}'])\n"
        "exit(0 if result == True or result is None else 1)\n"
    )

    execute_process(
        COMMAND ${LOCALPYTHON} ${EXTERNAL_ROOT}/install_FIAT.py
        RESULT_VARIABLE result
        ERROR_VARIABLE error
        OUTPUT_VARIABLE output
    )
    message(STATUS ${LOCALPYTHON})
    message(STATUS ${EXTERNAL_ROOT})
    message(STATUS ${ERROR_VARIABLE})
    message(STATUS ${OUTPUT_VARIABLE})
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
