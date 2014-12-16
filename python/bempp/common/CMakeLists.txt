set(headers __init__.pxd global_parameters.pxd)

install_python(FILES ${headers} 
               DESTINATION bempp/include/bempp/common)

add_python_module(bempp.common
    __init__.py ${headers} global_parameters.pyx
    TARGETNAME bempp.common
    CPP
    LIBRARIES libbempp
)

add_dependencies(cython-headers bempp.common)
