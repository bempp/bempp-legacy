set(makoes boundary_operator.mako.pxd boundary_operator.mako.pyx)
mako_files(${makoes}
    OUTPUT_FILES makoed
    DEPENDS "${PROJECT_SOURCE_DIR}/python/mako/space.py"
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/assembly"
    TARGETNAME bempp.assembly-mako
)

# Add empty __init__.pxd file so cython considers the include directory
set(init_pxd "${PYTHON_BINARY_DIR}/bempp/include/bempp/assembly/__init__.pxd")
file(WRITE "${init_pxd}")
list(APPEND makoed "${init_pxd}")

# copy other header file to include dir
configure_file(python.hpp
    "${PYTHON_BINARY_DIR}/bempp/include/bempp/assembly/python.hpp"
    COPYONLY
)

split_list(sources headers makoed ".*\\.pyx")
install_python(FILES ${headers} DESTINATION bempp/include/bempp/assembly)
add_dependencies(cython-headers bempp.assembly-mako)

add_python_module(bempp.assembly
    ${sources} __init__.py
    TARGETNAME bempp.assembly
    CPP
    LIBRARIES libbempp
)
add_dependencies(bempp.assembly cython-headers)
