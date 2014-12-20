mako_files(
    identity_operator.mako.pxd identity_operator.mako.pyx
    OUTPUT_FILES makoed
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/operators/sparse"
    TARGETNAME bempp.operators.sparse-mako
)
split_list(sources headers makoed ".*\\.pyx")

# Add empty __init__.pxd file so cython considers the include directory
set(init_pxd "${PYTHON_BINARY_DIR}/bempp/include/bempp/operators/sparse/__init__.pxd")
file(WRITE "${init_pxd}")

# copy other header file to include dir
list(APPEND headers "${init_pxd}")

install_python(FILES ${headers} DESTINATION bempp/include/bempp/operators/sparse)
add_dependencies(cython-headers bempp.operators.sparse-mako)

add_python_module(bempp.operators.sparse
    __init__.py ${sources}
    TARGETNAME bempp.operators.sparse
    CPP
    LIBRARIES libbempp
)
add_dependencies(bempp.operators.sparse cython-headers)
add_dependencies(cython-headers bempp.operators.sparse-mako)
