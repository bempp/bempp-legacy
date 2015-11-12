set(headers __init__.pxd
    )
            

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/core/operators/potential"
    TARGETNAME bempp.core.operators.potential-copied
)

install_python(FILES ${headers}
    DESTINATION bempp/include/bempp/core/operators/potential)
add_dependencies(cython-headers bempp.core.operators.potential-copied)

add_python_module(bempp.core.operators.potential __init__.py
                  laplace.pyx
                  modified_helmholtz.pyx
                  maxwell.pyx
                  CPP LIBRARIES libbempp)

add_dependencies(bempp.core.operators.potential bempp.core.assembly)
add_dependencies(bempp.core.operators.potential cython-headers)
add_dependencies(cython bempp.core.operators.potential)

