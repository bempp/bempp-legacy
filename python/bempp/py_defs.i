// Useful Python tools

%pythoncode %{

def checkType(dtype):
    dtypes={'float':'float64',
            'float32':'float32',
            'float64':'float64',
            'complex':'complex128',
            'complex64':'complex64',
            'complex128':'complex128'}
    if dtype in dtypes:
        return dtypes[dtype]
    else:
        raise ValueError('Data type does not exist.')

def promoteTypeToComplex(dtype):
    dtypes={'float':'complex128',
            'float32':'complex64',
            'float64':'complex128',
            'complex':'complex128',
            'complex64':'complex64',
            'complex128':'complex128'}
    if dtype in dtypes:
        return dtypes[dtype]
    else:
        raise ValueError('Data type does not exist.')

%}
