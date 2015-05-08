set(makoes __init__.mako.pxd)


mako_files(${makoes}
    OUTPUT_FILES makoed
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/file_interfaces/"
    TARGETNAME bempp.file_interfaces-mako
)


split_list(sources headers makoed ".*\\.pyx")


install_python(FILES ${headers} 
               DESTINATION bempp/include/bempp/file_interfaces)

add_python_module(bempp.file_interfaces
    __init__.py ${sources} 
    TARGETNAME bempp.file_interfaces
    CPP
    LIBRARIES libbempp
)

add_dependencies(cython-headers bempp.file_interfaces-mako)
add_dependencies(bempp.file_interfaces cython-headers)
