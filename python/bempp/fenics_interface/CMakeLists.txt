mako_files(
    _lagrange_coupling.mako.pxd _lagrange_coupling.mako.pyx
    py_lagrange_interface.mako.hpp
    OUTPUT_FILES makoed
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/fenics_interface"
    TARGETNAME bempp.fenics_interface-mako
    )

split_list(sources headers makoed ".*\\.pyx")

# Add empty __init__.pxd file so cython considers the include directory
set(init_pxd "${PYTHON_BINARY_DIR}/bempp/include/bempp/fenics_interface/__init__.pxd")
file(WRITE "${init_pxd}")

add_python_module(bempp.fenics_interface
    ${sources} __init__.py coupling.py p1_coupling.py fenics_operator.py
    TARGETNAME bempp.fenics_interface
    CPP
    LIBRARIES libbempp
)

add_dependencies(cython-headers bempp.fenics_interface-mako)
add_dependencies(bempp.fenics_interface cython-headers)
