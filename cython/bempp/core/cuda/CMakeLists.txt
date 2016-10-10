set(headers __init__.pxd)

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/core/cuda"
    TARGETNAME bempp.core.cuda-copied
)

install_python(FILES ${headers}
    DESTINATION bempp/include/bempp/core/cuda)
add_dependencies(cython-headers bempp.core.cuda-copied)

add_python_module(bempp.core.cuda __init__.py
                  test_cuda.pyx 
                  CPP LIBRARIES libbempp)

