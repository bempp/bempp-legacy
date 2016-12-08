find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()

# install and build paths for fake projects
set(PYTHON_BINARY_DIR "${PROJECT_BINARY_DIR}/python_binary"
    CACHE PATH "" FORCE)
set(PYTHON_PKG_DIR "${PROJECT_BINARY_DIR}/python_install"
    CACHE PATH "" FORCE)

find_package(CoherentPython)
include(PythonModule)

include(EnvironmentScript)
set(LOCAL_PYTHON_EXECUTABLE "@CMAKE_CURRENT_BINARY_DIR@/extension_tester.sh")
create_environment_script(
    EXECUTABLE "${PYTHON_EXECUTABLE}"
    PATH "${LOCAL_PYTHON_EXECUTABLE}"
    PYTHON
)
add_to_python_path("@EXTERNAL_ROOT@/python")
add_to_python_path("${PYTHON_BINARY_DIR}")

foreach(pathname module.c other.h other.cc)
    configure_file(@CMAKE_CURRENT_SOURCE_DIR@/${pathname} 
        "${CMAKE_CURRENT_SOURCE_DIR}/${pathname}"
        COPYONLY
    )
endforeach()
add_python_module("extension" module.c GLOB *.cc *.h)

foreach(target extension)
    if(NOT TARGET ${target})
        message(FATAL_ERROR "Target ${target} does not exist")
    endif()
    if(NOT TARGET ${target}-ext)
        message(FATAL_ERROR "Target ${target}-ext does not exist")
    endif()
endforeach()
