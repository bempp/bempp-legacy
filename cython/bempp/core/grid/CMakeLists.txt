set(files entity.pxd entity.pyx
           grid.pxd grid.pyx
           grid_view.pxd grid_view.pyx
           entity_iterator.pxd entity_iterator.pyx
           entity_pointer.pxd entity_pointer.pyx
           entity.pxd entity.pyx
           grid_factory.pyx
           geometry.pxd geometry.pyx
           codim_template.pxd
           index_set.pxd index_set.pyx
           py_sphere.hpp 
           py_entity_helper.hpp
           id_set.pxd id_set.pyx
           __init__.pxd)


split_list(sources headers files ".*\\.pyx")

copy_files_on_change(${headers}
           OUTPUT_FILES copied
           DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/core/grid"
           TARGETNAME bempp.core.grid-copied
           )

install_python(FILES ${headers} DESTINATION bempp/include/bempp/core/grid)
add_dependencies(cython-headers bempp.core.grid-copied)

add_python_module(bempp.core.grid __init__.py ${sources} CPP LIBRARIES libbempp)
add_dependencies(bempp.core.grid cython-headers)
add_dependencies(bempp.core.grid bempp.core.utils)
add_dependencies(cython bempp.core.grid)

