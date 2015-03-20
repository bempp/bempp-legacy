mako_files(helmholtz.mako.pyx maxwell.mako.pyx
    __init__.mako.pxd
    OUTPUT_FILES makoed
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/operators/far_field"
    TARGETNAME bempp.operators.far_field-mako
)

split_list(sources headers makoed ".*\\.pyx")

install_python(FILES ${headers} DESTINATION bempp/include/bempp/operators/far_field)

add_python_module(bempp.operators.far_field
    __init__.py ${sources}
    TARGETNAME bempp.operators.far_field
    CPP
    LIBRARIES libbempp
)
add_dependencies(cython-headers bempp.operators.far_field-mako)
add_dependencies(bempp.operators.far_field cython-headers)
add_dependencies(bempp.operators.far_field bempp.operators)

