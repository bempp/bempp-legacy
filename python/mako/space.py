dtypes = {
        'float32': 'float',
        'float64': 'double',
        'complex64': 'complex_float',
        'complex128': 'complex_double'
}
""" Possible dptypes for spaces and c equivalents """

class Class(object):
    def __init__(self, class_name):
        super(Class, self).__init__()
        self.class_name = class_name
    @property
    def type_name(self):
        return "%s[BASIS]" % self.class_name
    @property
    def header(self):
        f = lambda x: x if x.islower() else '_' + x.lower()
        name = ''.join([f(u) for u in self.class_name[1:]])
        return self.class_name[0].lower() + name

    scalar = property(lambda x: 'scalar' in x.class_name.lower())
    continuous = property(lambda x: 'discontinous' not in x.class_name.lower())
    constant = property(lambda x: 'constant' in x.class_name.lower())
    linear = property(lambda x: 'linear' not in x.class_name.lower())
    polynomial = property(lambda x: 'linear' not in x.class_name.lower())

# Contains all spaces
spaces = {}
def add(names, *args, **kwargs):
    """ Helper function for adding spaces """
    dictionary = spaces
    for name in names:
        if name not in dictionary:
            dictionary[name] = {}
        dictionary = dictionary[name]
    def adder(implementation='grid_only', **_kwargs):
        _kwargs['implementation'] = implementation
        dictionary[Class(*args, **kwargs)] = _kwargs
    return adder

# Now define individual implementations
add(
    ('scalar', 'constant', 'continuous'),
    'PiecewiseConstantScalarSpace'
)(
    doc='Space of piecewise constant scalar functions',
    tags=True
)

add(
    ('scalar', 'constant', 'continuous'),
    'PiecewiseConstantScalarSpaceBarycentric'
)(
    doc='Space of piecewise constant scalar functions',
    tags=True
)

add(
    ('scalar', 'constant', 'continuous'),
    'PiecewiseConstantDualGridScalarSpace'
)(
    doc='Space of piecewise constant scalar functions on the dual grid',
    tags=True,
)

add(
    ('scalar', 'linear', 'continuous'),
    'PiecewiseLinearContinuousScalarSpace'
)(
    doc='Space of continuous, piecewise linear scalar functions',
    tags=True
)


def flatten(spaces):
    from collections import Sequence
    for key, value in spaces.iteritems():
        if isinstance(key, Class):
            yield key, value
        else:
            for inners in flatten(value):
                yield inners

# Some sanity checks
# Class names are unique
class_names = [u[0].class_name for u in flatten(spaces)]
assert len(class_names) == len(set(class_names))
del class_names

