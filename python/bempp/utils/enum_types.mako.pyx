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

cdef TranspositionMode transposition_mode(string name):

    cdef TranspositionMode res

    if name==string(b'no_transpose'):
        res = no_transpose
    elif name==string(b'conjugate'):
        res = conjugate
    elif name==string(b'transpose'):
        res = transpose
    elif name==string(b'conjugate_transpose'):
        res = conjugate_transpose
    else:
        raise ValueError("Unsupported transposition mode")

    return res

cdef ConstructionMode construction_mode(string name):

    cdef ConstructionMode res

    if name==string(b'approximate'):
        res = approximate
    elif name==string(b'interpolate'):
        res = interpolate
    else:
        raise ValueError("Unsupported construction mode")

    return res


cdef GmshPostDataType gmsh_post_data_type(string name):

    cdef GmshPostDataType res

    if name==string(b'node'):
        res = node 
    elif name==string(b'element'):
        res = element
    elif name==string(b'element_node'):
        res = element_node
    else:
        raise ValueError("Unsupported gmsh type")

    return res
