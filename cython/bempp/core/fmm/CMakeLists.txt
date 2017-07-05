set(headers __init__.pxd
    octree.pxd)

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/core/fmm"
    TARGETNAME bempp.core.fmm-copied
    )

#install_python(FILES ${headers} DESTINATION bempp/include/bempp/core/fmm)
add_dependencies(cython-headers bempp.core.fmm-copied)

add_python_module(bempp.core.fmm __init__.py octree.pyx
    CPP LIBRARIES libbempp)
add_dependencies(bempp.core.fmm cython-headers)
add_dependencies(bempp.core.fmm bempp.core.grid)
add_dependencies(cython bempp.core.fmm)

