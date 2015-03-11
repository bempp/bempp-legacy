mako_files(
    identity_operator.mako.pxd identity_operator.mako.pyx maxwell_identity.mako.pxd
    maxwell_identity.mako.pyx
    OUTPUT_FILES makoed
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/operators/boundary/sparse"
    TARGETNAME bempp.operators.boundary.sparse-mako
)
split_list(sources headers makoed ".*\\.pyx")

# Add empty __init__.pxd file so cython considers the include directory
set(init_pxd "${PYTHON_BINARY_DIR}/bempp/include/bempp/operators/boundary/sparse/__init__.pxd")
file(WRITE "${init_pxd}")

# copy other header file to include dir
list(APPEND headers "${init_pxd}")

install_python(FILES ${headers} DESTINATION bempp/include/bempp/operators/boundary/sparse)
add_dependencies(cython-headers bempp.operators.boundary.sparse-mako)

add_python_module(bempp.operators.boundary.sparse
    __init__.py ${sources}
    TARGETNAME bempp.operators.boundary.sparse
    CPP
    LIBRARIES libbempp
)
add_dependencies(bempp.operators.boundary.sparse cython-headers)
add_dependencies(cython-headers bempp.operators.boundary.sparse-mako)
