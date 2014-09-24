set(makoes entity.mako.pxd entity.mako.pyx
           grid.mako.pxd grid.mako.pyx
           grid_view.mako.pxd grid_view.mako.pyx
           entity_iterator.mako.pxd entity_iterator.mako.pyx
           entity_pointer.mako.pxd entity_pointer.mako.pyx
           entity.mako.pxd entity.mako.pyx
           geometry.mako.pxd geometry.mako.pyx
           codim_template.mako.pxd)
mako_files(${makoes}
           OUTPUT_FILES makoed
           DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/grid"
           TARGETNAME bempp.grid-mako
           )


mako_files(__init__.mako.py
    OUTPUT_FILES init.py
    TARGETNAME bempp.grid-init-mako
)

# Add empty __init__.pxd file so cython considers the include directory
set(init_pxd "${PYTHON_BINARY_DIR}/bempp/include/bempp/grid/__init__.pxd")
file(WRITE "${init_pxd}")
list(APPEND makoed "${init_pxd}")

split_list(sources headers makoed ".*\\.pyx")
install_python(FILES ${headers} DESTINATION bempp/include/bempp/grid)
add_dependencies(cython-headers bempp.grid-mako)

add_python_module(bempp.grid ${init.py} ${sources} CPP LIBRARIES libbempp)
add_dependencies(bempp.grid cython-headers bempp.grid-init-mako)

# Add headers first
#add_python_module(bempp.grid 
#                  grid.pxd 
#                  grid_view.pxd 
#                  geometry.pxd
#                  entity_pointer.pxd
#                  entity_iterator.pxd
#                  codim_template.pxd
#                  entity.pxd TARGETNAME bempp.grid-headers)
#file(WRITE "${PYTHON_BINARY_DIR}/bempp/include/bempp/grid/__init__.pxd")

#add_python_module(bempp.grid 
#                  grid.pyx 
#                  grid_view.pyx 
#                  geometry.pyx
#                  entity_pointer.pyx
#                  entity_iterator.pyx
#                  entity.pyx
#                  __init__.py CPP LIBRARIES libbempp)
#
#add_dependencies(bempp.grid cython-headers)
#add_dependencies(cython-headers bempp.grid-headers)
