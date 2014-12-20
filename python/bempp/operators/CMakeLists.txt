mako_files( py_operators.mako.hpp 
    OUTPUT_FILES makoed
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/operators"
    TARGETNAME bempp.operators-mako
)

set(headers py_operators.mako.hpp)

# Add empty __init__.pxd file so cython considers the include directory
set(init_pxd "${PYTHON_BINARY_DIR}/bempp/include/bempp/operators/__init__.pxd")
file(WRITE "${init_pxd}")

# copy other header file to include dir
list(APPEND headers "${init_pxd}")

install_python(FILES ${headers} DESTINATION bempp/include/bempp/operators)
add_dependencies(cython-headers bempp.operators-mako)

add_python_module(bempp.operators
    __init__.py
    TARGETNAME bempp.operators
    CPP
    LIBRARIES libbempp
)
add_dependencies(bempp.operators cython-headers)
add_dependencies(cython-headers bempp.operators-mako)

add_subdirectory(sparse)

