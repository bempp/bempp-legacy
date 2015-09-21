set(headers __init__.pxd
    py_far_field_operators.hpp
    )
            

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp_ext/include/bempp_ext/operators/far_field"
    TARGETNAME bempp_ext.operators.far_field-copied
)

install_python(FILES ${headers}
    DESTINATION bempp_ext/include/bempp_ext/operators/far_field)
add_dependencies(cython-headers bempp_ext.operators.far_field-copied)

add_python_module(bempp_ext.operators.far_field __init__.py
                  helmholtz.pyx
                  maxwell.pyx
                  CPP LIBRARIES libbempp)

add_dependencies(bempp_ext.operators.far_field bempp_ext.assembly)
add_dependencies(bempp_ext.operators.far_field cython-headers)
