set(makoes __init__.mako.pxd)


mako_files(${makoes}
    OUTPUT_FILES makoed
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/external/"
    TARGETNAME bempp.external-mako
)

set(headers ${makoed})

install_python(FILES ${headers} 
               DESTINATION bempp/include/bempp/external)

add_python_module(bempp.external
    __init__.py ${sources} viewers.py
    TARGETNAME bempp.external
    CPP
    LIBRARIES libbempp
)

add_dependencies(cython-headers bempp.external-mako)
add_dependencies(bempp.external cython-headers)
