from py.test import mark
from bempp.assembly import Context


def test_scalar_space_creation():
    from bempp.grid import sample as sample_grid
    grid = sample_grid()

    context = Context()
    for dtype in ['float32', 'complex128']:
        context.basis_type = dtype
        constant = context.scalar_space(grid, order=0)
        assert constant.dtype == dtype


def test_dirichlet_tut_operators():
    from bempp.grid import sample as sample_grid
    grid = sample_grid()

    context = Context(basis_type='float32')
    constant = context.scalar_space(grid, order=0)
    linear = context.scalar_space(grid, order=1)

    double_layer = context.operators.laplace.boundary.double_layer(
        constant, linear, constant)
    single_layer = context.operators.laplace.boundary.single_layer(
        constant, linear, constant)
    identity = context.operators.local.identity(constant, linear, constant)

    assert double_layer.basis_type == context.basis_type
    assert single_layer.basis_type == context.basis_type
    assert identity.basis_type == context.basis_type

    assert single_layer.domain.is_compatible(constant)
    assert single_layer.range.is_compatible(linear)
    assert single_layer.dual_to_range.is_compatible(constant)


@mark.parametrize('result_type, scalar', [
    ('float32', 2.0),
    ('float64', 2.0),
    ('complex64', 0.5 + 2.0j),
    ('complex128', 0.5 + 2.0j),
    ('complex128', 0.5),
    # integer scalars should be promoted
    ('float64', 2),
    ('complex128', 2),
    ('complex128', 2 + 3j)
])
def test_boundary_op_operations(result_type, scalar):
    from py.test import raises
    from bempp.grid import sample as sample_grid
    grid = sample_grid()

    context = Context(result_type=result_type)
    constant = context.scalar_space(grid, order=0)
    linear = context.scalar_space(grid, order=1)

    double_layer = context.operators.laplace.boundary.double_layer(
        constant, linear, constant)
    single_layer = context.operators.laplace.boundary.single_layer(
        constant, linear, constant)

    sumop = double_layer + single_layer
    assert sumop != double_layer
    assert sumop != single_layer

    # Something about incompatible domains.
    # But that still checks that the operation is passed on to C++
    with raises(ValueError):
        double_layer * single_layer

    scalop = scalar * double_layer
    assert scalop != double_layer
    scalop = double_layer * scalar
    assert scalop != double_layer

    divop = double_layer / scalar
    assert divop != double_layer


@mark.parametrize('b0, b1, r0, r1', [
    ('float32', 'float32', 'float32', 'complex64'),
    ('float32', 'float64', 'float32', 'float64')
])
def test_incompatible_boundary_op_operations(b0, b1, r0, r1):
    from py.test import raises
    from bempp.grid import sample as sample_grid
    grid = sample_grid()

    c0 = Context(basis_type=b0, result_type=r0)
    c1 = Context(basis_type=b1, result_type=r1)

    op1 = c0.scalar_space(grid, order=0)
    op2 = c1.scalar_space(grid, order=0)
    with raises(TypeError):
        op1 + op2
    with raises(TypeError):
        op1 - op2
    with raises(TypeError):
        op1 * op2
