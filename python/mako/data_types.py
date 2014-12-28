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
    
def scalar_cython_type(name):

    return {
        'float': 'float',
        'float64': 'double',
        'complex_float': 'float complex',
        'complex_double': 'double complex'
        }.get(name,name)

def real_cython_type(name):

    return {
        'float': 'float',
        'float64': 'double',
        'complex_float': 'float',
        'complex_double': 'double'
        }.get(name,name)
