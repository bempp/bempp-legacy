set(copy_files py_types.hpp py_utils.hpp
    parameter_list.pxd parameter_list.pyx
    __init__.pxd eigen.pxd eigen.pyx shared_ptr.pxd
    unique_ptr.pxd enum_types.pxd enum_types.pyx
    signal_slot_interface.pxd)

copy_files_on_change(${copy_files}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/utils/"
    TARGETNAME bempp.utils-copied
)


split_list(sources headers copied ".*\\.pyx")


install_python(FILES ${headers} 
               DESTINATION bempp/include/bempp/utils)

add_python_module(bempp.utils
    __init__.py data_types.py which.py matrix_operator.py
    blocked_linear_operator.py byte_conversion.py ${sources} 
    TARGETNAME bempp.utils
    CPP
    LIBRARIES libbempp
)

add_dependencies(cython-headers bempp.utils-copied)
add_dependencies(bempp.utils cython-headers)
