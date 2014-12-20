<%!
from data_types import dtypes
from space import spaces
def get_key(description):
    return (
        description['order'],
        description['continuous'],
        description['barycentric'],
        description['dual'],
        description.get('extra', None),
    )
%>
__all__ = ['Space', 'scalar']
from .space import Space
from . import space as _cwrappers


_scalar_spaces = {
    # (order, continuous, barycentric, dual, extra): implementation
% for class_name, description in spaces.items():
%   if description['scalar']:
    ${get_key(description)}: _cwrappers.${class_name},
%   endif
% endfor
}
""" All scalar spaces with definitions """


def scalar_class(barycentric=False, continuous=True, order=None,
           constant=None, linear=None, dual=False, extra=None):
<%def name="tags_params()" filter="trim">
        barycentric : bool
            Whether the space is barycentric. Defaults to False.

        order: 0|1|>=2
            Order of the functions. Order 0 (1) is equivalent to the constant
            (linear) keyword.

        constant: bool
            Constant piecewise functions.
        linear: bool
            Linear piecewise functions.

        dual: bool
            If True, then the functions are applied on the dual grid.

        continuous: bool
            If True, the space consists of continous functions.
</%def>\
<%def name="all_spaces()" filter="trim">
<%
    from operator import itemgetter
    def row(order, continuous, dual, barycentric, **kwargs):
        return "{order: <10}  {continuous: ^10}  "\
                "{dual: ^9}  {barycentric: ^11}".format(
            order=order,
            continuous='yes' if continuous else 'no',
            dual='yes' if dual else 'no',
            barycentric='yes' if barycentric else 'no',
        )
%>\
        The default space is the space of continuous piecewise-constant
        functions on the direct grid. The following set of spaces are
        available:

        ==========  ==========  =========  ===========
        order       continuous  dual grid  barycentric
        ==========  ==========  =========  ===========
%   for description in sorted(spaces.values(), key=itemgetter('order')):
        ${row(**description)}
%   endfor
        ==========  ==========  =========  ===========
</%def>\
    """ Discriminates between spaces according to input

        Parameters
        ----------

        ${tags_params()}

        ${all_spaces()}
    """
    # input sanity check and sanitization
    def nargs(*args):
        sum([int(bool(u)) for u in args])
    kind = None
    # default call
    if order is None and constant is None and linear is None:
        kind = 'constant'
    # non-default calls
    elif (order is None and constant) or order == 0 or order == 'constant':
        kind = 'constant'
    if (order is None and linear) or order == 1 or order == 'linear':
        if kind is not None:
            raise TypeError("Input requests both %s and linear" % kind)
        kind = 'linear'
    if kind is None:
        kind = 'linear' if order is None else 'polynomial'

    # Now figure out which to use
    key = (kind, continuous, barycentric, dual, extra)
    if key not in _scalar_spaces:
        raise ValueError("No space exists for this set of input combinations")
    return _scalar_spaces[key]


def scalar(grid, dtype, barycentric=False, continuous=True, order=None,
           constant=None, linear=None, dual=False, extra=None, **kwargs):
    """ Factory for creating scalar spaces

        Parameters
        ----------

        grid : Grid
            Grid over which the functions are applied

        dtype : ${'|'.join(dtypes.keys())}
            Precision and kind of the functions.

        ${tags_params()}

        Other keyword parameters are passed on when constructing the spaces per
        se.

        ${all_spaces()}
    """
    # passes order on to constructor
    if order is not None and order > 1:
        kwargs['order'] = order
    return scalar_class(
        barycentric=barycentric, continuous=continuous, order=order,
        constant=constant, linear=linear, dual=dual, extra=extra
    )(grid, dtype, **kwargs)
