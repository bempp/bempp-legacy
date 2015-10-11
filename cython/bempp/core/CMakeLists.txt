set(headers __init__.pxd)

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/core"
    TARGETNAME bempp.core-copied
)

install_python(FILES __init__.pxd DESTINATION bempp/include/bempp/core)
add_dependencies(cython-headers bempp.core-copied)

add_python_module(bempp.core __init__.py CPP LIBRARIES libbempp)
add_dependencies(bempp.core cython-headers)

add_subdirectory(utils)
add_subdirectory(common)
add_subdirectory(grid)
add_subdirectory(space)
add_subdirectory(assembly)
add_subdirectory(operators)
add_subdirectory(hmat)
