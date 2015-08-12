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
           DESTINATION "${PYTHON_BINARY_DIR}/bempp_ext/include/bempp_ext/grid"
           TARGETNAME bempp_ext.grid-copied
           )

install_python(FILES ${headers} DESTINATION bempp/include/bempp/grid)
add_dependencies(cython-headers bempp_ext.grid-copied)

add_python_module(bempp_ext.grid __init__.py ${sources} CPP LIBRARIES libbempp)
add_dependencies(bempp_ext.grid cython-headers)
