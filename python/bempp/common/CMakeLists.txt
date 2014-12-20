set(makoes __init__.mako.pxd global_parameters.mako.pxd global_parameters.mako.pyx)

mako_files(${makoes}
           OUTPUT_FILES makoed
           DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/common"
           TARGETNAME bempp.common-mako)

split_list(sources headers makoed ".*\\.pyx")

install_python(FILES ${headers} 
               DESTINATION bempp/include/bempp/common)
add_dependencies(cython-headers bempp.common-mako)

add_python_module(bempp.common
    __init__.py ${sources}
    TARGETNAME bempp.common
    CPP
    LIBRARIES libbempp
)

add_dependencies(bempp.common cython-headers)
