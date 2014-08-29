from py.test import mark, fixture
from bempp.space.space import PiecewiseConstantScalarSpace, \
    PiecewiseLinearContinuousScalarSpace


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
