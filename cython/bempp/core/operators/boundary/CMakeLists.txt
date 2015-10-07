set(headers __init__.pxd
    support_operators.hpp
    )
            

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/core/operators/boundary"
    TARGETNAME bempp.core.operators.boundary-copied
)

install_python(FILES ${headers}
    DESTINATION bempp/include/bempp/core/operators/boundary)
add_dependencies(cython-headers bempp.core.operators.boundary-copied)

add_python_module(bempp.core.operators.boundary __init__.py
                  laplace.pyx
                  sparse.pyx
                  modified_helmholtz.pyx
                  maxwell.pyx
                  CPP LIBRARIES libbempp)

add_dependencies(bempp.core.operators.boundary cython-headers)
