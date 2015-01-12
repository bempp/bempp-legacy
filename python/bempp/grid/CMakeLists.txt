set(makoes entity.mako.pxd entity.mako.pyx
           grid.mako.pxd grid.mako.pyx
           grid_view.mako.pxd grid_view.mako.pyx
           entity_iterator.mako.pxd entity_iterator.mako.pyx
           entity_pointer.mako.pxd entity_pointer.mako.pyx
           entity.mako.pxd entity.mako.pyx
           geometry.mako.pxd geometry.mako.pyx
           codim_template.mako.pxd
           py_sphere.mako.hpp)
mako_files(${makoes}
           OUTPUT_FILES makoed
           DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/grid"
           TARGETNAME bempp.grid-mako
           )


# Add empty __init__.pxd file so cython considers the include directory
set(init_pxd "${PYTHON_BINARY_DIR}/bempp/include/bempp/grid/__init__.pxd")
file(WRITE "${init_pxd}")
list(APPEND makoed "${init_pxd}")

split_list(sources headers makoed ".*\\.pyx")
install_python(FILES ${headers} DESTINATION bempp/include/bempp/grid)
add_dependencies(cython-headers bempp.grid-mako)

add_python_module(bempp.grid __init__.py ${sources} CPP LIBRARIES libbempp)
add_dependencies(bempp.grid cython-headers)
