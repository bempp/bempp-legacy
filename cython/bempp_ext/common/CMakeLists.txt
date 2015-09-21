set(headers __init__.pxd)

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp_ext/include/bempp_ext/common"
    TARGETNAME bempp_ext.common-copied
)

install_python(FILES __init__.pxd DESTINATION bempp_ext/include/bempp_ext/common)
add_dependencies(cython-headers bempp_ext.common-copied)

add_python_module(bempp_ext.common __init__.py global_parameters.pyx CPP LIBRARIES libbempp)
add_dependencies(bempp_ext.common cython-headers)
