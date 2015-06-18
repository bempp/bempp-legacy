add_python_module(bempp.shapes
    shapes.py __init__.py 
    TARGETNAME bempp.shapes
    CPP
    LIBRARIES libbempp
)

add_dependencies(bempp.shapes cython-headers)
