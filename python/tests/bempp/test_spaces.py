from py.test import mark, fixture
from bempp.space.space import PiecewiseConstantScalarSpace, \
    PiecewiseLinearContinuousScalarSpace,                   \
    PiecewiseLinearDiscontinuousScalarSpace,                \
    PiecewisePolynomialContinuousScalarSpace,               \
    PiecewiseConstantDualGridScalarSpace


@fixture
def grid():
    from os.path import join, exists
    from bempp.config import paths
    from bempp.grid import Grid
    filename = join(paths.meshes, "sphere-h-0.4.msh")
    if not exists(filename):
        raise IOError("Mesh %s does not exist" % filename)
    return Grid(topology="triangular", filename=filename)


@mark.parametrize("TestClass, dtype", [
    (PiecewiseConstantScalarSpace, 'float32'),
    (PiecewiseConstantScalarSpace, 'float64'),
    (PiecewiseConstantScalarSpace, 'complex64'),
    (PiecewiseConstantScalarSpace, 'complex128'),
    (PiecewiseLinearContinuousScalarSpace, 'float64')
])
def test_instantiation(grid, TestClass, dtype):
    space = TestClass(grid, dtype)
    assert space.dtype == dtype
    assert space.grid is not None
    assert space.grid == grid


def test_grid_is_readonly(grid):
    from py.test import raises
    space = PiecewiseLinearContinuousScalarSpace(grid, 'complex128')
    with raises(AttributeError):
        space.grid = grid


@mark.parametrize("TestClass, kwargs", [
    (PiecewiseConstantScalarSpace, {}),
    (PiecewiseConstantScalarSpace, {'constant': True}),
    (PiecewiseConstantScalarSpace, {'order': 0}),
    (PiecewiseConstantScalarSpace, {'order': 'constant'}),
    (PiecewiseLinearDiscontinuousScalarSpace,
        {'linear': True, 'continuous': False}),
    (PiecewiseLinearDiscontinuousScalarSpace,
        {'order': 1, 'continuous': False}),
    (PiecewiseLinearDiscontinuousScalarSpace,
        {'order': 'linear', 'continuous': False}),
    (PiecewisePolynomialContinuousScalarSpace, {'order': 2}),
    (PiecewisePolynomialContinuousScalarSpace, {'order': 3})
])
def test_space_selection(TestClass, kwargs):
    from bempp.space import scalar_class

    actual = scalar_class(**kwargs)
    assert actual is TestClass


def test_polynomial_order(grid):
    from bempp.space import scalar

    for order in [2, 3]:
        actual = scalar(grid, 'float32', order=order)
        assert isinstance(actual, PiecewisePolynomialContinuousScalarSpace)
        assert actual.order == order


def test_space_is_compatible(grid):
    constant64 = PiecewiseConstantScalarSpace(grid, 'complex64')
    constant128 = PiecewiseConstantScalarSpace(grid, 'complex128')

    assert not constant64.is_compatible(constant128)
    assert not constant128.is_compatible(constant64)
    assert constant64.is_compatible(constant64)

    linear64 = PiecewiseLinearContinuousScalarSpace(grid, 'complex64')
    linear64b = PiecewiseLinearContinuousScalarSpace(grid, 'complex64')
    assert not constant64.is_compatible(linear64)
    assert linear64b.is_compatible(linear64)
    assert linear64.is_compatible(linear64b)

    dual = PiecewiseConstantDualGridScalarSpace(grid, 'complex64')
    assert not constant64.is_compatible(dual)
