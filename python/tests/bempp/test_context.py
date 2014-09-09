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

    double_layer = context.operators.laplace_3d.double_layer(
        constant, linear, constant)
    single_layer = context.operators.laplace_3d.single_layer(
        constant, linear, constant)
    identity = context.operators.identity(constant, linear, constant)

    assert double_layer.basis_type == context.basis_type
    assert single_layer.basis_type == context.basis_type
    assert identity.basis_type == context.basis_type
