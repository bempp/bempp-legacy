set(headers py_types.hpp py_utils.hpp
    parameter_list.pxd parameter_list.pyx
    __init__.pxd eigen.pxd eigen.pyx shared_ptr.pxd
    unique_ptr.pxd enum_types.pxd enum_types.pyx
    signal_slot_interface.pxd complex.pxd)

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp_ext/include/bempp_ext/utils/"
    TARGETNAME bempp_ext.utils-copied
)


install_python(FILES ${headers} 
               DESTINATION bempp_ext/include/bempp_ext/utils)

add_python_module(bempp_ext.utils
    parameter_list.pyx eigen.pyx enum_types.pyx __init__.py
    TARGETNAME bempp_ext.utils
    CPP
    LIBRARIES libbempp
)

add_dependencies(cython-headers bempp_ext.utils-copied)
add_dependencies(bempp_ext.utils cython-headers)
