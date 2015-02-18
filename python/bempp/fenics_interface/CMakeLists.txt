set(sources __init__.py coupling.py)
add_python_module(bempp.fenics_interface
    ${sources}
    TARGETNAME bempp.fenics_interface
    CPP
    LIBRARIES libbempp
)

add_dependencies(bempp.fenics_interface cython-headers)
