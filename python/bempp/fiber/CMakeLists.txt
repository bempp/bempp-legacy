add_python_module("bempp.fiber"
    quadrature_options.pyx quadrature_options.pxd __init__.py
    CPP
    LIBRARIES bempp
    TARGETNAME python_bempp
)
