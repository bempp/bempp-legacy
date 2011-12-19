include_directories(${PROJECT_SOURCE_DIR}/lib)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs)
find_package(PythonInterp)
find_package(Numpy)
include_directories(${PYTHON_INCLUDE_PATH})
# ".." added to remove the "/numpy" ending from ${PYTHON_NUMPY_INCLUDE_DIR}
include_directories(${PYTHON_NUMPY_INCLUDE_DIR}/..)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Remove "-modern" if you need compatibility with Python 2.1.x and older
set(CMAKE_SWIG_FLAGS "-modern")

set_source_files_properties(core.i PROPERTIES CPLUSPLUS ON)
set(SWIG_MODULE_core_EXTRA_DEPS core.i armadillo2numpy.i auto_ptr.i numpy.i 
    docstrings.i
    grid/entity.i grid/entity_iterator.i grid/entity_pointer.i 
    grid/geometry.i grid/geometry_type.i
    grid/grid.i grid/grid_factory.i grid/grid_parameters_converter.hpp grid/grid_view.i
    grid/id_set.i grid/index_set.i grid/vtk_writer.i
    grid/entity_pointer_docstrings.i grid/entity_iterator_docstrings.i
    grid/grid_docstrings.i grid/grid_view_docstrings.i 
    grid/id_set_docstrings.i grid/index_set_docstrings.i 
    grid/geometry_docstrings.i grid/geometry_type_docstrings.i 
    grid/grid_factory_docstrings.i grid/vtk_writer_docstrings.i
)
swig_add_module(core python core.i)
swig_link_libraries(core ${PYTHON_LIBRARIES} bempp)

### Move __init__.py to the build directory ###################################

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

###############################################################################