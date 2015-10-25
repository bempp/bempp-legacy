set(headers __init__.pxd)

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/core/common"
    TARGETNAME bempp.core.common-copied
)

install_python(FILES __init__.pxd DESTINATION bempp/include/bempp/core/common)
add_dependencies(cython-headers bempp.core.common-copied)

add_python_module(bempp.core.common __init__.py global_parameters.pyx CPP LIBRARIES libbempp)
add_dependencies(bempp.core.common cython-headers)
add_dependencies(cython bempp.core.common)
