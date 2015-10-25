set(headers assembler.pxd
            abstract_boundary_operator.pxd 
            discrete_boundary_operator.pxd
            discrete_operator_conversion.hpp
            function_projector.hpp
            __init__.pxd
            )

copy_files_on_change(${headers}
    OUTPUT_FILES copied
    DESTINATION "${PYTHON_BINARY_DIR}/bempp/include/bempp/core/assembly"
    TARGETNAME bempp.core.assembly-copied
)

install_python(FILES ${headers}
    DESTINATION bempp/include/bempp/core/assembly)
add_dependencies(cython-headers bempp.core.assembly-copied)

add_python_module(bempp.core.assembly __init__.py
                  assembler.pyx
                  abstract_boundary_operator.pyx
                  discrete_boundary_operator.pyx 
                  function_projector.pyx
                  CPP LIBRARIES libbempp)

add_dependencies(bempp.core.assembly cython-headers)
add_dependencies(bempp.core.assembly bempp.core.grid)
add_dependencies(bempp.core.assembly bempp.core.utils)
add_dependencies(bempp.core.assembly bempp.core.space)
add_dependencies(cython bempp.core.assembly)
