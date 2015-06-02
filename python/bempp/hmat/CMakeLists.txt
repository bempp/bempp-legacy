mako_files(
    hmatrix.mako.pxd hmatrix.mako.pyx block_cluster_tree.mako.pxd
    block_cluster_tree.mako.pyx hmatrix_data.mako.pxd 
    hmatrix_data.mako.pyx py_hmat_support.mako.hpp
    OUTPUT_FILES makoed
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/hmat"
    TARGETNAME bempp.hmat-mako
)
split_list(sources headers makoed ".*\\.pyx")

# Add empty __init__.pxd file so cython considers the include directory
set(init_pxd "${PYTHON_BINARY_DIR}/bempp/include/bempp/hmat/__init__.pxd")
file(WRITE "${init_pxd}")

# copy other header file to include dir
list(APPEND headers "${init_pxd}")

install_python(FILES ${headers} DESTINATION bempp/include/bempp/hmat)
add_dependencies(cython-headers bempp.hmat-mako)

add_python_module(bempp.hmat
    ${sources} __init__.py 
    TARGETNAME bempp.hmat
    CPP
    LIBRARIES libbempp
)
add_dependencies(bempp.hmat cython-headers)
add_dependencies(cython-headers bempp.hmat-mako)
