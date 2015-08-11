copy_files_on_change(space.pxd space.pyx __init__.pxd
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/space"
    TARGETNAME bempp.space-copied
)

split_list(sources headers copied ".*\\.pyx")
install_python(FILES ${headers} DESTINATION bempp/include/bempp/space)
add_dependencies(cython-headers bempp.space-copied)

add_python_module(bempp.space __init__.py ${sources} CPP LIBRARIES libbempp)
add_dependencies(bempp.space cython-headers)
