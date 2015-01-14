
mako_files(laplace.mako.pxd laplace.mako.pyx modified_helmholtz.mako.pxd 
    modified_helmholtz.mako.pyx helmholtz.mako.pyx
    OUTPUT_FILES makoed
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/operators/boundary"
    TARGETNAME bempp.operators.boundary-mako
)

split_list(sources headers makoed ".*\\.pyx")

# Add empty __init__.pxd file so cython considers the include directory
set(init_pxd "${PYTHON_BINARY_DIR}/bempp/include/bempp/operators/boundary/__init__.pxd")
file(WRITE "${init_pxd}")

# copy other header file to include dir
list(APPEND headers "${init_pxd}")

install_python(FILES ${headers} DESTINATION bempp/include/bempp/operators/boundary)

add_python_module(bempp.operators.boundary
    __init__.py ${sources}
    TARGETNAME bempp.operators.boundary
    CPP
    LIBRARIES libbempp
)
add_dependencies(bempp.operators.boundary cython-headers)
add_dependencies(bempp.operators.boundary bempp.operators)
add_dependencies(cython-headers bempp.operators.boundary-mako)
add_subdirectory(sparse)

