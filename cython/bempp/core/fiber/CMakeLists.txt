set(headers __init__.pxd shapeset.pxd
            arrays.pxd)

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/core/fiber"
    TARGETNAME bempp.core.fiber-copied
)

install_python(FILES ${headers}
    DESTINATION bempp/include/bempp/core/fiber)
add_dependencies(cython-headers bempp.core.fiber-copied)

add_python_module(bempp.core.fiber __init__.py
                  integrators.pyx arrays.pyx shapeset.pyx
                  CPP LIBRARIES libbempp)

add_dependencies(bempp.core.fiber cython-headers)
add_dependencies(bempp.core.fiber bempp.core.utils)
add_dependencies(cython bempp.core.fiber)
