set(headers __init__.pxd)
            

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/core/operators"
    TARGETNAME bempp.core.operators-copied
)

install_python(FILES ${headers}
    DESTINATION bempp/include/bempp/core/operators)
add_dependencies(cython-headers bempp.core.operators-copied)

add_python_module(bempp.core.operators __init__.py
                  CPP LIBRARIES libbempp)

add_dependencies(bempp.core.operators cython-headers)
add_dependencies(bempp.core.operators bempp.core.assembly)
add_dependencies(cython bempp.core.operators)

add_subdirectory(boundary)
add_subdirectory(potential)
add_subdirectory(far_field)
