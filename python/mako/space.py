dtypes = {
    'float32': 'float',
    'float64': 'double',
    'complex64': 'complex_float',
    'complex128': 'complex_double'
}
""" Possible dptypes for spaces and c equivalents """

compatible_dtypes = {  # name: (is real, precision)
    'float32': ('float32', 'complex64'),
    'float64': ('float64', 'complex128'),
    'complex64': ('complex64',),
    'complex128': ('complex128',),
}
""" Compatibility between basis and result types """


# Describes available spaces and their wrapper implementation.
# Most of the characteristics (barycentric, dual...) are guessed later on.
# However, these characteristics can be specified here. The characteristics
# will be used in the actual space factories by the users.
# The default implementation is 'grid_only', eg a single constructor that takes
# a grid as its only argument.
spaces = {
    'PiecewiseConstantScalarSpace': {
        'doc': 'Space of piecewise constant scalar functions',
    },
    'PiecewiseConstantScalarSpaceBarycentric': {
        'doc': 'Space of piecewise constant scalar functions',
    },
    'PiecewiseConstantDualGridScalarSpace': {
        'doc': 'Space of piecewise constant scalar functions on the dual grid',
    },
    'PiecewiseLinearContinuousScalarSpace': {
        'doc': 'Space of continuous, piecewise linear scalar functions',
    },
    'PiecewiseLinearDiscontinuousScalarSpace': {
        'doc':
        'Space of piecewise linear, possibly discontinuous, scalar functions',
    },
    'PiecewiseLinearDiscontinuousScalarSpaceBarycentric': {
        'doc':
        'Space of piecewise constant functions define on the dual grid',
        'dual': True  # According top the doxygen tag...
    },
    'PiecewisePolynomialContinuousScalarSpace': {
        'doc':
        'Space of continuous, piecewise polynomial scalar functions',
        'implementation': 'polynomial'
    }
}


# Guess implementation from class name
# Mostly means guessing wether space operates on the direct or dual grid,
# whether functions are continuous, whether they are constant, linear, or
# polynomial, etc. These facts are used later on in the actual space factory.
for key, description in spaces.iteritems():
    if 'implementation' not in description:
        description['implementation'] = 'grid_only'

    if 'header' not in description:
        f = lambda x: x if x.islower() else '_' + x.lower()
        description['header'] = 'bempp/space/%s.hpp' \
            % (key[0].lower() + ''.join([f(u) for u in key[1:]]))

    if 'scalar' not in description:
        description['scalar'] = 'scalar' in key.lower()

    if description['scalar'] and 'order' not in description:
        if 'constant' in key.lower():
            description['order'] = 'constant'
        elif 'linear' in key.lower():
            description['order'] = 'linear'
        else:
            description['order'] = 'polynomial'

    if 'continuous' not in description:
        description['continuous'] = 'discontinuous' not in key.lower()

    if 'dual' not in description:
        description['dual'] = 'dual' in key.lower()

    if 'barycentric' not in description:
        description['barycentric'] = 'barycentric' in key.lower()
