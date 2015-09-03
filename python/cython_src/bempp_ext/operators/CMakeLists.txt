set(headers __init__.pxd)
            

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp_ext/include/bempp_ext/operators"
    TARGETNAME bempp_ext.operators-copied
)

install_python(FILES ${headers}
    DESTINATION bempp_ext/include/bempp_ext/operators)
add_dependencies(cython-headers bempp_ext.operators-copied)

add_python_module(bempp_ext.operators __init__.py
                  CPP LIBRARIES libbempp)

add_dependencies(bempp_ext.operators cython-headers)

add_subdirectory(boundary)
add_subdirectory(potential)
