find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()

include(PythonPackageLookup)
include(EnvironmentScript)
include(MakoFiles)

set(LOCAL_PYTHON_EXECUTABLE "${PROJECT_BINARY_DIR}/localpython.sh")
create_environment_script(
    EXECUTABLE "${PYTHON_EXECUTABLE}"
    PATH "${LOCAL_PYTHON_EXECUTABLE}"
    PYTHON
)
add_to_python_path("${EXTERNAL_ROOT}/python")

lookup_python_package(mako REQUIRED)
find_program(mako_SCRIPT mako-render HINT "${EXTERNAL_ROOT}/python")
if(NOT mako_SCRIPT)
    message(FATAL_ERROR "Could not find mako-render script.")
endif()

file(WRITE "${CMAKE_CURRENT_SOURCE_DIR}/__init__.mako.py"
    "import other\n"
    "i = 0\n"
    "% for a in ['hello', 'world']:\n"
    "assert '\${a}' == 'hello world'.split()[\${loop.index}]\n"
    "i += 1\n"
    "% endfor\n"
    "assert i == 2\n"
    "assert other.i == 8\n"
)
file(WRITE "${CMAKE_CURRENT_SOURCE_DIR}/other.mako.py"
    "i = 5\n"
    "% for a in ['hello', 'despicable', 'world']:\n"
    "assert '\${a}' == 'hello despicable world'.split()[\${loop.index}]\n"
    "i += 1\n"
    "% endfor\n"
    "assert i == 8\n"
)

set(destination "${CMAKE_CURRENT_BINARY_DIR}/python_package/makoed")
mako_files(GLOB *.mako.py
    DESTINATION "${destination}"
    OUTPUT_FILES output
)
# output must be used somewhere. Otherwise it is not built.
# In practice, the output will be used in some library or something.
# Otherwise, a new target with the output files need to be declared as below.
# Doing add_dependencies on a prior target seems to fail, however.
add_custom_target(makoed ALL DEPENDS ${output})

list(LENGTH output i)
if(NOT i EQUAL 2)
    message(FATAL_ERROR "Expected 2 output files. Got ${output}.")
endif()

foreach(filename other.py __init__.py)
    set(filename "${destination}/${filename}")
    list(FIND output "${filename}" found)
    if(found LESS 0)
        message(FATAL_ERROR "${filename} not in output files")
    endif()
endforeach()

install(FILES ${output} DESTINATION
    "${CMAKE_CURRENT_BINARY_DIR}/install")
