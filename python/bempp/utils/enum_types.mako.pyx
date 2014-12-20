from libcpp.string cimport string

cdef SymmetryMode symmetry_mode(string name):

    cdef SymmetryMode res

    if name==string(b'no_symmetry'):
        res = no_symmetry
    elif name==string(b'symmetric'):
        res = symmetric
    elif name==string(b'hermitian'):
        res = hermitian
    elif name==string(b'auto_symmetry'):
        res = auto_symmetry
    else:
        raise ValueError("Unsupported symmetry mode")

    return res

