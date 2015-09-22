set(headers space.pxd __init__.pxd local_evaluator.hpp)

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/core/space"
    TARGETNAME bempp.core.space-copied
)

install_python(FILES ${headers} DESTINATION bempp/include/bempp/core/space)
add_dependencies(cython-headers bempp.core.space-copied)

add_python_module(bempp.core.space __init__.py space.pyx CPP LIBRARIES libbempp)
add_dependencies(bempp.core.space cython-headers)
