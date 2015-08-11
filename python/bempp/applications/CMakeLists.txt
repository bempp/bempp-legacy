add_python_module(bempp.applications
    __init__.py helmholtz.py
    TARGETNAME bempp.applications
    CPP
    LIBRARIES libbempp
)

add_dependencies(bempp.applications cython-headers)
