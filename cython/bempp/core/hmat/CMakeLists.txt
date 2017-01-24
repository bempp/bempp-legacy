set(headers __init__.pxd
            block_cluster_tree.pxd
            py_hmat_support.hpp
            hmatrix_data.pxd)

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/core/hmat"
    TARGETNAME bempp.core.hmat-copied
)

#install_python(FILES ${headers} DESTINATION bempp/include/bempp/core/hmat)
add_dependencies(cython-headers bempp.core.hmat-copied)

add_python_module(bempp.core.hmat __init__.py block_cluster_tree.pyx 
                                  hmatrix_data.pyx hmatrix_interface.pyx
                                  CPP LIBRARIES libbempp)
add_dependencies(bempp.core.hmat cython-headers)
add_dependencies(bempp.core.hmat bempp.core.assembly)
add_dependencies(cython bempp.core.hmat)

