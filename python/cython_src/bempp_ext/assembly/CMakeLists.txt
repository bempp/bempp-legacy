set(headers abstract_boundary_operator.pxd 
            discrete_boundary_operator.pxd
            __init__.pxd
            )

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp_ext/include/bempp_ext/assembly"
    TARGETNAME bempp_ext.assembly-copied
)

install_python(FILES ${headers}
    DESTINATION bempp_ext/include/bempp_ext/assembly)
add_dependencies(cython-headers bempp_ext.assembly-copied)

add_python_module(bempp_ext.assembly __init__.py
                  abstract_boundary_operator.pyx
                  discrete_boundary_operator.pyx 
                  CPP LIBRARIES libbempp)

add_dependencies(bempp_ext.assembly cython-headers)
