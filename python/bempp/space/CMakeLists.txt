set(makoes space.mako.pxd space.mako.pyx types.mako.h __init__.mako.py)
mako_files(${makoes}
    OUTPUT_FILES makoed
    DEPENDS "${PROJECT_SOURCE_DIR}/python/mako/space.py"
)

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/__init__.pxd")
list(APPEND makod "${CMAKE_CURRENT_BINARY_DIR}/__init__.pxd")
split_list(headers sources makoed ".*\\.pxd")

add_custom_target(bempp.space-headers DEPENDS ${headers})
add_python_module(bempp.space ${headers} types.h TARGETNAME bempp.space-headers)
add_dependencies(cython-headers bempp.space-headers)

add_python_module(bempp.space ${sources} CPP LIBRARIES libbempp)
add_dependencies(bempp.space bempp.space-headers)

# CMake does not always pick up dependencies from one add_custom_command to
# another. So we add the dependency explicitly
add_custom_target(bempp.space-mako DEPENDS ${makoed})
add_dependencies(bempp.space-headers bempp.space-mako)
