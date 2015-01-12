mako_files(
    laplace.pxd laplace.pyx
    OUTPUT_FILES makoed
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/operators/laplace"
    TARGETNAME bempp.operators.laplace-mako
)
split_list(sources headers makoed ".*\\.pyx")

# Add empty __init__.pxd file so cython considers the include directory
set(init_pxd "${PYTHON_BINARY_DIR}/bempp/include/bempp/operators/laplace/__init__.pxd")
file(WRITE "${init_pxd}")

# copy other header file to include dir
list(APPEND headers "${init_pxd}")

install_python(FILES ${headers} DESTINATION bempp/include/bempp/operators/laplace)
add_dependencies(cython-headers bempp.operators.laplace-mako)

add_python_module(bempp.operators.laplace
    __init__.py ${sources}
    TARGETNAME bempp.operators.laplace
    CPP
    LIBRARIES libbempp
)
add_dependencies(bempp.operators.laplace cython-headers)
add_dependencies(cython-headers bempp.operators.laplace-mako)

