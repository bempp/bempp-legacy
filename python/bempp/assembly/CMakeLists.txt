mako_files(boundary_operator.mako.pyx python.mako.hpp context.mako.pyx
    OUTPUT_FILES makoed
    DEPENDS "${PROJECT_SOURCE_DIR}/python/mako/space.py"
        "${PROJECT_SOURCE_DIR}/python/mako/bempp_operators.py"
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/assembly"
    TARGETNAME bempp.assembly-mako
)
split_list(sources headers makoed ".*\\.pyx")

# Add empty __init__.pxd file so cython considers the include directory
set(init_pxd "${PYTHON_BINARY_DIR}/bempp/include/bempp/assembly/__init__.pxd")
file(WRITE "${init_pxd}")

# copy other header file to include dir
list(APPEND headers "${init_pxd}")
foreach(filename boundary_operator.pxd)
    configure_file(${filename}
        "${PYTHON_BINARY_DIR}/bempp/include/bempp/assembly/${filename}"
        COPYONLY
    )
    list(APPEND headers ${filename})
endforeach()

install_python(FILES ${headers} DESTINATION bempp/include/bempp/assembly)
add_dependencies(cython-headers bempp.assembly-mako)

add_python_module(bempp.assembly
    ${sources} __init__.py
    TARGETNAME bempp.assembly
    CPP
    LIBRARIES libbempp
)
add_dependencies(bempp.assembly cython-headers)
add_dependencies(cython-headers bempp.assembly-mako)
