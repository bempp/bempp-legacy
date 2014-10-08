set(makoes space.mako.pxd space.mako.pyx types.mako.h variant.mako.hpp)
mako_files(${makoes}
    OUTPUT_FILES makoed
    DEPENDS "${PROJECT_SOURCE_DIR}/python/mako/space.py"
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/space"
    TARGETNAME bempp.space-mako
)

mako_files(__init__.mako.py
    OUTPUT_FILES init.py
    DEPENDS "${PROJECT_SOURCE_DIR}/python/mako/space.py"
    TARGETNAME bempp.space-init-mako
)

# Add empty __init__.pxd file so cython considers the include directory
set(init_pxd "${PYTHON_BINARY_DIR}/bempp/include/bempp/space/__init__.pxd")
file(WRITE "${init_pxd}")
list(APPEND makoed "${init_pxd}")

split_list(sources headers makoed ".*\\.pyx")
install_python(FILES ${headers} DESTINATION bempp/include/bempp/space)
add_dependencies(cython-headers bempp.space-mako)

add_python_module(bempp.space ${init.py} ${sources} CPP LIBRARIES libbempp)
add_dependencies(bempp.space cython-headers bempp.space-init-mako)
