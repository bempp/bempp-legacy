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

set_source_files_properties(bempp.i PROPERTIES CPLUSPLUS ON)
set(SWIG_MODULE_bempp_EXTRA_DEPS bempp.i armadillo2numpy.i auto_ptr.i numpy.i 
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
swig_add_module(bempp python bempp.i)
swig_link_libraries(bempp ${PYTHON_LIBRARIES} bempp)
