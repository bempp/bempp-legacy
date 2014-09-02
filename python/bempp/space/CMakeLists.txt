set(makoes space.mako.pxd space.mako.pyx types.mako.h __init__.mako.py)
mako_files(${makoes}
    OUTPUT_FILES makoed
    DEPENDS "${PROJECT_SOURCE_DIR}/python/mako/space.py"
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/space"
    TARGETNAME bempp.space-mako
)

# Add empty __init__.pxd file so cython considers the directory
set(init_pxd "${PYTHON_BINARY_DIR}/bempp/include/bempp/space/__init__.pxd")
file(WRITE "${init_pxd}")
list(APPEND makoed "${init_pxd}")

split_list(sources headers makoed ".*\\.pyx")
install_python(FILES ${headers} DESTINATION bempp/include/bempp/space)
add_dependencies(cython-headers bempp.space-mako)

add_python_module(bempp.space ${sources} CPP LIBRARIES libbempp)
add_dependencies(bempp.space cython-headers)
