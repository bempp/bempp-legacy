add_python_module("bempp.grid"
    grid.pyx grid.pxd __init__.py
    EXCLUDE config.in.py
    CPP
    LIBRARIES bempp
    TARGETNAME python_bempp
)
