# include swig functionality
include(${SWIG_USE_FILE})
# Remove "-modern" if you need compatibility with Python 2.1.x and older
set(CMAKE_SWIG_FLAGS "-modern")

include_directories(${PYTHON_INCLUDE_DIR} ${NUMPY_INCLUDE_DIRS})

set_source_files_properties(core.i PROPERTIES CPLUSPLUS ON)
# Use the following to debug swig
# set_source_files_properties(core.i PROPERTIES SWIG_FLAGS "-debug-classes")
file(GLOB_RECURSE SWIG_SOURCES *.i)
set(SWIG_MODULE_core_EXTRA_DEPS ${SWIG_SOURCES})
swig_add_module(core python core.i)
swig_link_libraries(core ${PYTHON_LIBRARY} bempp teuchoscore teuchosparameterlist)

# Postprocess docstrings
add_custom_command(
    TARGET ${SWIG_MODULE_core_REAL_NAME}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} copy -E
        ${CMAKE_CURRENT_BINARY_DIR}/core.py
        ${CMAKE_CURRENT_BINARY_DIR}/core.py.orig
    COMMAND ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/postprocess_docstrings.py
        ${CMAKE_CURRENT_BINARY_DIR}/core.py
)
configure_file(
	${CMAKE_SOURCE_DIR}/python/bempp/pyfragments.swg
    ${CMAKE_BINARY_DIR}/python/bempp/pyfragments.swg)

### Move __init__.py to the build directory ###################################
if (NOT (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR}))
    #set(PYTHON_SOURCES "__init__.py")
    file(GLOB_RECURSE PYTHON_SOURCES RELATIVE
        ${CMAKE_SOURCE_DIR}/python/bempp
        ${CMAKE_SOURCE_DIR}/python/bempp/*.py
    )

    foreach(py_file ${PYTHON_SOURCES})
        add_custom_command(
            OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${py_file}
            COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/${py_file}
                ${CMAKE_CURRENT_BINARY_DIR}/${py_file}
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${py_file}
            COMMENT "Copying file ${CMAKE_CURRENT_SOURCE_DIR}/${py_file}"
            VERBATIM
        )

	    set(PYTHON_DEPS "${PYTHON_DEPS};${CMAKE_CURRENT_BINARY_DIR}/${py_file}")
    endforeach()

    add_custom_target(copy_python_bempp ALL DEPENDS ${PYTHON_DEPS})
endif ()

###############################################################################

file(GLOB SWIG_SWG_FILES *.swg)
file(GLOB PYTHON_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.py)
list(REMOVE_ITEM PYTHON_FILES "determineNumpyPath.py" "postprocess_docstrings.py")
install_python(TARGETS _core DESTINATION bempp)
install_python(FILES ${PYTHON_FILES}
    "${CMAKE_CURRENT_BINARY_DIR}/core.py" DESTINATION bempp)
install_python(FILES ${SWIG_SWG_FILES} ${SWIG_SOURCES} core.i DESTINATION bempp/include)
