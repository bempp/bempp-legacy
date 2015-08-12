set(headers space.pxd __init__.pxd)

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp_ext/include/bempp_ext/space"
    TARGETNAME bempp_ext.space-copied
)

install_python(FILES space.pxd __init__.pxd DESTINATION bempp_ext/include/bempp_ext/space)
add_dependencies(cython-headers bempp_ext.space-copied)

add_python_module(bempp_ext.space __init__.py space.pyx CPP LIBRARIES libbempp)
add_dependencies(bempp_ext.space cython-headers)
