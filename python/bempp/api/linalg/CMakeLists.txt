install_python(FILES __init__.pxd 
               DESTINATION bempp/include/bempp/linalg)

add_python_module(bempp.api.linalg
    __init__.py iterative_solvers.py direct_solvers.py blocked.py
    TARGETNAME bempp.api.linalg
    CPP
    LIBRARIES libbempp
)

add_dependencies(bempp.api.linalg cython-headers)
