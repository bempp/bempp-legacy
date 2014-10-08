set(headers quadrature_options.pxd range_accuracy_options.pxd
    accuracy_options.pxd)
add_python_module(bempp.fiber ${headers} TARGETNAME bempp.fiber-headers)
file(WRITE "${PYTHON_BINARY_DIR}/bempp/include/bempp/fiber/__init__.pxd")

add_python_module(bempp.fiber
    __init__.py quadrature_options.pyx range_accuracy_options.pyx
    accuracy_options.pyx
    CPP
    LIBRARIES libbempp
)
add_dependencies(bempp.fiber cython-headers)
add_dependencies(cython-headers bempp.fiber-headers)
