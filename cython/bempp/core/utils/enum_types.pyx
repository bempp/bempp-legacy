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


cdef HMatBlockType hmat_block_type(string name):

    cdef HMatBlockType res

    if name==string(b'dense'):
        res = dense 
    elif name==string(b'low_rank_ab'):
        res = low_rank_ab 
    else:
        raise ValueError("Unsupported block type")

    return res

cdef TranspositionMode compute_transpose_mode(
        TranspositionMode current_mode,
        TranspositionMode input_mode):

    if current_mode==no_transpose:
        return input_mode
    if current_mode==transpose:
        if input_mode==no_transpose:
            return transpose
        if input_mode==transpose:
            return no_transpose
        if input_mode==conjugate:
            return conjugate_transpose
        if input_mode==conjugate_transpose:
            return transpose
    if current_mode==conjugate:
        if input_mode==no_transpose:
            return conjugate
        if input_mode==transpose:
            return conjugate_transpose
        if input_mode==conjugate:
            return no_transpose
        if input_mode==conjugate_transpose:
            return transpose
    if current_mode==conjugate_transpose:
        if input_mode==no_transpose:
            return conjugate_transpose
        if input_mode==transpose:
            return conjugate
        if input_mode==conjugate:
            return transpose
        if input_mode==conjugate_transpose:
            return no_transpose

       



