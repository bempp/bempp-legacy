set(headers __init__.pxd
    py_potential_operators.hpp
    )
            

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp_ext/include/bempp_ext/operators/potential"
    TARGETNAME bempp_ext.operators.potential-copied
)

install_python(FILES ${headers}
    DESTINATION bempp_ext/include/bempp_ext/operators/potential)
add_dependencies(cython-headers bempp_ext.operators.potential-copied)

add_python_module(bempp_ext.operators.potential __init__.py
                  laplace.pyx
                  modified_helmholtz.pyx
                  maxwell.pyx
                  CPP LIBRARIES libbempp)

add_dependencies(bempp_ext.operators.potential bempp_ext.assembly)
add_dependencies(bempp_ext.operators.potential cython-headers)
