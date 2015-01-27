import numpy as np


dtypes = {
    'float32': 'float',
    'float64': 'double',
    'complex64': 'complex_float',
    'complex128': 'complex_double'
}
""" Possible dptypes for spaces and c equivalents """

compatible_dtypes = {  # name: (tuple of compatible dtypes)
    'float32': ('float32', 'complex64'),
    'float64': ('float64', 'complex128'),
    'complex64': ('complex64',),
    'complex128': ('complex128',),
}
""" Compatibility between basis and result types """
def ctypes(name):
    """ Valid c++ type if name is a cython or python type 

        Meant only to work with types from dtypes and compatible_dtypes
    """
    return {
        'float32': 'float',
        'float64': 'double',
        'complex64': 'std::complex<float>',
	'complex64': 'std::complex<float>',
	'complex_float': 'std::complex<float>',
        'complex_double': 'std::complex<double>'
    }.get(name, name)
    
def check_type(name,default='float64'):
    """ Try to convert input into a numpy.dtype object """


    value = str(name) if name is not None else default
    if value not in ['float32', 'float64', 'complex64', 'complex128']:
        raise ValueError("Incorrect type (%s)" % value)

    return np.dtype(value)

def combined_type(dtype1,dtype2):
    """ Return a type that is compatible with dtype1 and dtype2 """

    import numpy as np
    d1 = check_type(dtype1)
    d2 = check_type(dtype2)

    return (d1.type(1)*d2.type(1)).dtype



