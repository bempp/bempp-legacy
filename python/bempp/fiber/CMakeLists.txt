add_python_module("bempp.fiber"
    __init__.py
    quadrature_options.pyx quadrature_options.pxd
    range_accuracy_options.pyx range_accuracy_options.pxd
    accuracy_options.pyx accuracy_options.pxd
    CPP
    LIBRARIES bempp
    TARGETNAME python_bempp
)
