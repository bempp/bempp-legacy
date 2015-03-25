install_python(FILES __init__.pxd 
               DESTINATION bempp/include/bempp/linalg)

add_python_module(bempp.linalg
    __init__.py iterative_solvers.py direct_solvers.py
    TARGETNAME bempp.linalg
    CPP
    LIBRARIES libbempp
)

add_dependencies(bempp.linalg cython-headers)
