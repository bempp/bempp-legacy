<%!
from space import spaces, Class

def module_names(spaces, name=''):
    for key, value in spaces.iteritems():
        if not isinstance(key, Class):
            new_name = key if name == '' else name + "." + key
            yield new_name
            for name_ in module_names(value, new_name):
                yield name_
%>
__all__ = ['Space', 'scalar']
from .space import Space

def scalar(grid, dtype, barycentric=False, continuous=True, order=None,
        constant=False, linear=False, polynomial=False, **kwargs):
    """ Creates scalar space """
    # input sanity check and sanitization
    if order is not None:
        order = int(order)
        if order < 0:
            raise ValueError("order should positive or null")
        condition = (order == 0 and (constant is None or constant)) \
            or (order == 1 and (linear is None or linear)) \
            or (polynomial is None or polynomial))
        if not condition:
            raise ValueError("Inconsistent arguments")
        if order > 2:
            order = 2
    elif sum([constant, linear, polynomial]) > 1:
        raise ValueError("Inconsistent arguments")
    elif sum([constant, linear, polynomial]) = 0 or constant:
        order = 0
    elif linear:
        linear = 1
    else:
        order = 2


    # Now figure out which to use
    possible = ${possibles}
    key = str(order) + str(barycentric) + str(continuous)
    if key not in possibles:
        raise ValueError("No space exists for this set of input combinations")
    return possibles[key](grid, dtype, *args, **kwargs)
