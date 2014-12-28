set(makoes space.mako.pxd space.mako.pyx py_space_variants.mako.hpp __init__.mako.pxd __init__.mako.pyx)
mako_files(${makoes}
    OUTPUT_FILES makoed
    DEPENDS "${PROJECT_SOURCE_DIR}/python/mako/space.py"
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/space"
    TARGETNAME bempp.space-mako
)

split_list(sources headers makoed ".*\\.pyx")
install_python(FILES ${headers} DESTINATION bempp/include/bempp/space)
add_dependencies(cython-headers bempp.space-mako)

add_python_module(bempp.space ${sources} CPP LIBRARIES libbempp)
add_dependencies(bempp.space cython-headers)
