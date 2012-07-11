include_directories(${PROJECT_SOURCE_DIR}/lib)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs)
find_package(PythonInterp)
find_package(Numpy)
include_directories(${PYTHON_INCLUDE_PATH})
# ".." added to remove the "/numpy" ending from ${PYTHON_NUMPY_INCLUDE_DIR}
include_directories(${PYTHON_NUMPY_INCLUDE_DIR}/..)

include_directories ("${CMAKE_SOURCE_DIR}/lib")
include_directories(${CMAKE_BINARY_DIR}/include/common)
include_directories(${ARMADILLO_INCLUDE_DIR})
include_directories(${BOOST_INCLUDE_DIR})
include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
include_directories(${CMAKE_INSTALL_PREFIX}/bempp/contrib/dune/include)


# Remove "-modern" if you need compatibility with Python 2.1.x and older
set(CMAKE_SWIG_FLAGS "-modern")

set_source_files_properties(core.i PROPERTIES CPLUSPLUS ON)
file(GLOB_RECURSE SWIG_SOURCES *.i)
set(SWIG_MODULE_core_EXTRA_DEPS ${SWIG_SOURCES})
swig_add_module(core python core.i)
swig_link_libraries(core ${PYTHON_LIBRARIES} bempp)

### Move __init__.py to the build directory ###################################

if (NOT (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR}))
    set(PYTHON_SOURCES "__init__.py")
    foreach (py_file ${PYTHON_SOURCES})
    add_custom_command(
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${py_file}
	COMMAND ${CMAKE_COMMAND} copy -E ${CMAKE_CURRENT_SOURCE_DIR}/${py_file} ${CMAKE_CURRENT_BINARY_DIR}/${py_file}
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${py_file}
	# TODO: remove the project's root directory from the printed path.
	COMMENT "Copying file ${CMAKE_CURRENT_SOURCE_DIR}/${py_file}"
	VERBATIM
    )
    endforeach ()

    foreach(py_file ${PYTHON_SOURCES})
	set(PYTHON_DEPS "${PYTHON_DEPS};${CMAKE_CURRENT_BINARY_DIR}/${py_file}")
    endforeach()
    add_custom_target(copy_python_bempp ALL DEPENDS ${PYTHON_DEPS})
endif ()

###############################################################################
