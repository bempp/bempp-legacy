# Finds python interpreter and library.
# Tries to make sure one fits the other
if(PYTHON_INCLUDE_DIR AND PYTHON_LIBRARIES)
    return()
endif()
unset(required)
if(CoherentPython_FIND_REQUIRED)
    set(required REQUIRED)
endif()
unset(quietly)
if(CoherentPython_FIND_QUIETLY)
    set(quietly QUIET)
endif()
unset(version)
if(CoherentPython_FIND_VERSION)
    set(version ${CoherentPython_FIND_VERSION})
endif()
unset(exact)
if(CoherentPython_FIND_VERSION_EXACT)
    set(exact EXACT)
endif()

# Finds python interpreter first
find_package(PythonInterp ${version} ${exact} ${required} ${quietly})
if(NOT PYTHONINTERP_FOUND)
    return()
endif()

# Finds prefix from python executable itself, then adds
# relevant paths to relevant cmake variables so that
# find_package(PythonLibs) actually works
if(NOT CMAKE_CROSS_COMPILING)
    include(CallPython)
    # Figures out prefix from sys.prefix variable
    # Add relevant path to relevant variables, including MACOSX framework.
    macro(add_to_framework_paths)
        call_python(PYTHON_INTERP_PREFIX "import sys; print(sys.prefix)")
        if(DEFINED PYTHON_INTERP_PREFIX)
            list(INSERT CMAKE_PREFIX_PATH 0 ${PYTHON_INTERP_PREFIX})
            if("${PYTHON_INTERP_PREFIX}" MATCHES ".*/Frameworks/.*")
                string(REGEX REPLACE
                    "(.*/Frameworks)/.*" "\\1"
                    framework_path
                    "${PYTHON_INTERP_PREFIX}"
                )
                list(INSERT CMAKE_FRAMEWORK_PATH 0 "${framework_path}")
            endif()
        endif()
    endmacro()
    # Figures out paths from distutils' sysconfig module and sys.prefix
    # Started off for canopy and virtualenv
    macro(paths_from_distutils_sysconfig)
        call_python(PYTHON_INTERP_PREFIX "import sys; print(sys.prefix)")
        if(DEFINED PYTHON_INTERP_PREFIX)
            FILE(TO_CMAKE_PATH ${PYTHON_INTERP_PREFIX} PYTHON_INTERP_PREFIX)
            if(EXISTS "${PYTHON_INTERP_PREFIX}/libs")
                # Conda stores its .lib files in ${PYTHON_INTERP_PREFIX/libs
                # on windows. (Python35.lib etc needed for linking)
                # The dlls should be on path.
                list(INSERT CMAKE_LIBRARY_PATH 0 "${PYTHON_INTERP_PREFIX}/libs")
            endif()
        endif()
        call_python(python_include
        "from distutils.sysconfig import get_python_inc"
            "print(get_python_inc())"
        )
        if(DEFINED python_include AND EXISTS "${python_include}")
            FILE(TO_CMAKE_PATH ${python_include} python_include)
            set(PYTHON_INCLUDE_DIR "${python_include}")
        endif()
        # And tries adding sysconfig.get_python_lib output
        # Good for canopy, but not virtualenv.
        call_python(python_lib
        "from distutils.sysconfig import get_python_lib"
            "print(get_python_lib())"
        )
        if(DEFINED python_lib)
            FILE(TO_CMAKE_PATH ${python_lib} python_lib)
            if(EXISTS "${python_lib}")
                list(INSERT CMAKE_LIBRARY_PATH 0 "${python_lib}")
            endif()
            if(EXISTS "${python_lib}/../config")
                # CASA seems to put stuff here
                list(INSERT CMAKE_LIBRARY_PATH 0 "${python_lib}/../config")
            endif()
        endif()
        # And tries get_config_var('LIBDIR').
        # Good for virtualenv but not canopy.
        call_python(python_lib
        "from distutils.sysconfig import get_config_var"
        "print(get_config_var('LIBDIR'))"
        )
        if(DEFINED python_lib AND EXISTS "${python_lib}")
            list(INSERT CMAKE_LIBRARY_PATH 0 "${python_lib}")
        endif()
    endmacro()

    add_to_framework_paths()
    paths_from_distutils_sysconfig()
endif()

set(version_string "${PYTHON_VERSION_MAJOR}")
set(version_string "${version_string}.${PYTHON_VERSION_MINOR}")
set(version_string "${version_string}.${PYTHON_VERSION_PATCH}")
find_package(PythonLibs ${version_string} EXACT ${required} ${quiet})
if(NOT PYTHON_INCLUDE_DIRS)
    message(FATAL_ERROR "Could not find python include directory")
endif()
