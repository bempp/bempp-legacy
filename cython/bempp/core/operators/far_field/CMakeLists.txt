set(headers __init__.pxd
    )
            

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/core/operators/far_field"
    TARGETNAME bempp.core.operators.far_field-copied
)

install_python(FILES ${headers}
    DESTINATION bempp/include/bempp/core/operators/far_field)
add_dependencies(cython-headers bempp.core.operators.far_field-copied)

add_python_module(bempp.core.operators.far_field __init__.py
                  helmholtz.pyx
                  maxwell.pyx
                  CPP LIBRARIES libbempp)

add_dependencies(bempp.core.operators.far_field bempp.core.assembly)
add_dependencies(bempp.core.operators.far_field cython-headers)
add_dependencies(cython bempp.core.operators.far_field)

