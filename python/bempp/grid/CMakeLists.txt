# Add headers first
add_custom_target(bempp.grid-headers DEPENDS grid.pxd)
add_python_module(bempp.grid grid.pxd TARGETNAME bempp.grid-headers)
file(WRITE "${PYTHON_BINARY_DIR}/bempp/include/bempp/grid/__init__.pxd")

add_python_module(bempp.grid grid.pyx __init__.py CPP LIBRARIES libbempp)
add_dependencies(bempp.grid cython-headers)
add_dependencies(cython-headers bempp.grid-headers)
