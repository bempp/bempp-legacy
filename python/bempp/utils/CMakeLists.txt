file(COPY utils.hpp DESTINATION ${PYTHON_BINARY_DIR}/bempp/include/bempp/utils)

set(headers __init__.pxd utils.hpp armadillo.pxd parameter_list.pxd shared_ptr.pxd unique_ptr.pxd)

install_python(FILES ${headers} 
               DESTINATION bempp/include/bempp/utils)

add_python_module(bempp.utils
    __init__.py byte_conversion.py ${headers} parameter_list.pyx
    TARGETNAME bempp.utils
    CPP
    LIBRARIES libbempp
)

add_dependencies(cython-headers bempp.utils)
