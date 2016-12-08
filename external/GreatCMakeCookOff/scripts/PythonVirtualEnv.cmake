# Creates a python virtual environment in the build directory
# See https://github.com/UCL/GreatCMakeCookOff/wiki for information
find_package(PythonInterp)
include(Utilities)
include(PythonPackage)
include(EnvironmentScript)

if(NOT VIRTUALENV_WAS_CREATED)
    function(_create_virtualenv_from_exec call)
        execute_process(COMMAND
            ${call} ${VIRTUALENV_DIRECTORY}
            RESULT_VARIABLE RESULT
            ERROR_VARIABLE ERROR
            OUTPUT_VARIABLE OUTPUT
        )
        if(NOT "${RESULT}" STREQUAL "0")
            message(STATUS "${RESULT}")
            message(STATUS "${OUTPUT}")
            message(STATUS "${ERROR}")
            message(FATAL_ERROR "Could not create virtual environment")
        endif()
    endfunction()
    function(_create_virtualenv_from_package PACKAGE)
        execute_process(COMMAND
            ${PYTHON_EXECUTABLE} -m ${PACKAGE} ${VIRTUALENV_DIRECTORY}
            RESULT_VARIABLE RESULT
            ERROR_VARIABLE ERROR
            OUTPUT_VARIABLE OUTPUT
        )
        if(NOT "${RESULT}" STREQUAL "0")
            message(STATUS "${RESULT}")
            message(STATUS "${OUTPUT}")
            message(STATUS "${ERROR}")
            message(FATAL_ERROR "Could not create virtual environment")
        endif()
    endfunction()

    # First check if we have venv in the same directory as python.
    # Canopy makes things more difficult.
    get_filename_component(python_bin "${PYTHON_EXECUTABLE}" PATH)
    find_program(venv_EXECUTABLE venv PATHS "${python_bin}" NO_DEFAULT_PATH)
    find_python_package(venv QUIET)
    find_python_package(virtualenv QUIET)

    set(VIRTUALENV_DIRECTORY "${CMAKE_BINARY_DIR}/external/virtualenv"
        CACHE INTERNAL "Path to virtualenv" )
    if(venv_EXECUTABLE)
      _create_virtualenv_from_exec("${venv_EXECUTABLE}")
    elseif(venv_FOUND)
      _create_virtualenv_from_package("venv")
    elseif(virtualenv_FOUND)
      _create_virtualenv_from_package("virtualenv")
    else()
      message(FATAL_ERROR "Could find neither venv nor virtualenv")
    endif()

    set(_LOCAL_PYTHON_EXECUTABLE
        "${VIRTUALENV_DIRECTORY}/bin/python"
        CACHE PATH
        "Path to python executable in build virtual environment"
    )
    # Adds a bash script to call python with all the right paths set
    # Makes it easy to debug and test
    set(LOCAL_PYTHON_EXECUTABLE
        "${PROJECT_BINARY_DIR}/localpython.sh"
        CACHE PATH "Path to proxy for executing python script in build dir"
    )
    create_environment_script(
        EXECUTABLE "${_LOCAL_PYTHON_EXECUTABLE}"
        PATH "${PROJECT_BINARY_DIR}/localpython.sh"
    )

    # Add current python paths to a path.pth file
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -c
           "from sys import path;print(';'.join([u for u in path if len(u)]))"
        ${PROJECT_BINARY_DIR}/CMakeFiles/pypaths.py
        OUTPUT_VARIABLE output
    )
    foreach(pypath ${output})
        file(APPEND "${PROJECT_BINARY_DIR}/paths/system.pth" "${pypath}\n")
    endforeach()
    # Makes sure the path.pth file is picked up
    execute_process(COMMAND
        ${_LOCAL_PYTHON_EXECUTABLE} -c "import site; print(site.__file__)"
        OUTPUT_VARIABLE sitedir
        ERROR_VARIABLE error
        RESULT_VARIABLE result
    )
    if("${result}" STREQUAL "0")
        string(STRIP sitedir "${sitedir}")
        get_filename_component(sitedir "${sitedir}" PATH)
        file(WRITE "${sitedir}/sitecustomize.py"
            "from site import addsitedir\n"
            "addsitedir('${PROJECT_BINARY_DIR}/paths')\n"
        )
    else()
        message("error: ${error}")
        message("out: ${sitedir}")
        message(FATAL_ERROR "script failed")
    endif()

    set(VIRTUALENV_WAS_CREATED TRUE CACHE INTERNAL "Virtualenv has been created")
endif()

function(add_package_to_virtualenv PACKAGE)
    find_python_package(${PACKAGE} LOCAL)
    if(NOT ${PACKAGE}_FOUND)
        # First check if we have venv in the same directory as python.
        # Canopy makes things more difficult.
        get_filename_component(python_bin "${_LOCAL_PYTHON_EXECUTABLE}" PATH)
        find_program(local_pip_EXECUTABLE pip PATHS "${python_bin}" NO_DEFAULT_PATH)
        if(local_pip_EXECUTABLE)
            execute_process(
                COMMAND
                    ${local_pip_EXECUTABLE} install --upgrade ${PACKAGE}
                RESULT_VARIABLE result
                OUTPUT_VARIABLE output
                ERROR_VARIABLE error
            )
        else()
            execute_process(
                COMMAND ${_LOCAL_PYTHON_EXECUTABLE} -m pip install --upgrade ${PACKAGE}
                RESULT_VARIABLE result
                OUTPUT_VARIABLE output
                ERROR_VARIABLE error
            )
        endif()
        if("${result}" STREQUAL "0")
            find_python_package(${PACKAGE} LOCAL)
        else()
            message("${error}")
            message("${output}")
            message(FATAL_ERROR "Could not install ${PACKAGE} -- ${result}")
        endif()
    endif()
endfunction()
