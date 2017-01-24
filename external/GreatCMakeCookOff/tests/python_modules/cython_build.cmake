find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()

# install and build paths for fake projects
set(PYTHON_BINARY_DIR "${PROJECT_BINARY_DIR}/python_binary"
    CACHE PATH "" FORCE)
set(PYTHON_PKG_DIR "${PROJECT_BINARY_DIR}/python_install"
    CACHE PATH "" FORCE)

find_package(CoherentPython)
include(PythonModule)
include(PythonPackageLookup)
include(EnvironmentScript)

set(LOCAL_PYTHON_EXECUTABLE "@CMAKE_CURRENT_BINARY_DIR@/cython_tester.sh")
create_environment_script(
    EXECUTABLE "${PYTHON_EXECUTABLE}"
    PATH "${LOCAL_PYTHON_EXECUTABLE}"
    PYTHON
)
add_to_python_path("@EXTERNAL_ROOT@/python")
add_to_python_path("${PYTHON_BINARY_DIR}")

lookup_python_package(pytest REQUIRED PATH "@EXTERNAL_ROOT@/python")
lookup_python_package(cython REQUIRED PATH "@EXTERNAL_ROOT@/python")
get_filename_component(directory "${PYTHON_EXECUTABLE}" PATH)
find_program(cython_EXECUTABLE cython HINTS "${directory}")

foreach(pathname structure.pyx structure.pxd structure.h structure.c)
    configure_file("@CMAKE_CURRENT_SOURCE_DIR@/${pathname}"
        "${CMAKE_CURRENT_SOURCE_DIR}/${pathname}"
        COPYONLY
    )
endforeach()

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
add_library(pystructure SHARED structure.c)
add_python_module("extension" GLOB *.pyx *.pxd LIBRARIES pystructure FAKE_INIT)
