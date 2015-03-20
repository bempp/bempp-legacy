mako_files(
    py_operators.mako.hpp
    OUTPUT_FILES makoed
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/operators"
    TARGETNAME bempp.operators-mako
    )

# Add empty __init__.pxd file so cython considers the include directory
set(init_pxd "${PYTHON_BINARY_DIR}/bempp/include/bempp/operators/__init__.pxd")
file(WRITE "${init_pxd}")

install_python(FILES "${init_pxd}" DESTINATION bempp/include/bempp/operators)

add_python_module(bempp.operators
    __init__.py 
    TARGETNAME bempp.operators
    CPP
    LIBRARIES libbempp
    )
add_dependencies(cython-headers bempp.operators-mako)
add_dependencies(bempp.operators cython-headers)
add_subdirectory(boundary)
add_subdirectory(potential)
add_subdirectory(far_field)
