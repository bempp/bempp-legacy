mako_files(laplace.mako.pyx modified_helmholtz.mako.pyx helmholtz.mako.pyx maxwell.mako.pyx
    __init__.mako.pxd
    OUTPUT_FILES makoed
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/operators/potential"
    TARGETNAME bempp.operators.potential-mako
)

split_list(sources headers makoed ".*\\.pyx")

install_python(FILES ${headers} DESTINATION bempp/include/bempp/operators/potential)

add_python_module(bempp.operators.potential
    __init__.py ${sources}
    TARGETNAME bempp.operators.potential
    CPP
    LIBRARIES libbempp
)
add_dependencies(cython-headers bempp.operators.potential-mako)
add_dependencies(bempp.operators.potential cython-headers)
add_dependencies(bempp.operators.potential bempp.operators)

