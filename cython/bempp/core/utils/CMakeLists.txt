set(headers py_types.hpp py_utils.hpp
    parameter_list.pxd parameter_list.pyx
    __init__.pxd eigen.pxd eigen.pyx shared_ptr.pxd
    unique_ptr.pxd enum_types.pxd enum_types.pyx
    signal_slot_interface.pxd complex.pxd)

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/core/utils/"
    TARGETNAME bempp.core.utils-copied
)


install_python(FILES ${headers} 
               DESTINATION bempp/include/bempp/core/utils)

# STARTLINE adds a line a the beginning of the generated cython files
# It makes sure that Python.h is included *after* property tree
# This is necessary on mac (or maybe clang++) to avoid unpleasant conflicts
# between ctypes.h and cctypes.
add_python_module(bempp.core.utils
    parameter_list.pyx eigen.pyx enum_types.pyx __init__.py
    TARGETNAME bempp.core.utils
    CPP
    LIBRARIES libbempp
    STARTLINE "#include <boost/property_tree/ptree.hpp>"
)

add_dependencies(cython-headers bempp.core.utils-copied)
add_dependencies(bempp.core.utils cython-headers)
