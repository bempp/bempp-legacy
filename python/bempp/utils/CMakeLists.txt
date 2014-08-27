add_python_module(bempp.utils
    __init__.pxd armadillo.pxd utils.h shared_ptr.pxd)
add_dependencies(cython-headers bempp.utils)
