import dolfin as _dolfin

def boundary_grid_from_fenics_mesh(fenics_mesh):
    """
    Return a boundary grid from a FEniCS Mesh.
    """
    from bempp import grid_from_element_data

    bm = _dolfin.BoundaryMesh(fenics_mesh, "exterior", False)
    bm_coords = bm.coordinates()
    bm_cells  = bm.cells()
    bempp_boundary_grid = grid_from_element_data(bm_coords.transpose(),bm_cells.transpose())
    return bempp_boundary_grid

def generic_pn_trace(fenics_space):

    import dolfin
    from .coupling import fenics_space_info
    from bempp import function_space, grid_from_element_data, GridFunction, global_parameters
    import numpy as np

    family,degree = fenics_space_info(fenics_space)
    if not (family=='Lagrange'):
        raise ValueError("fenics_space different from Lagrange space not yet tested!")

    mesh = fenics_space.mesh()

    bm = dolfin.BoundaryMesh(mesh,"exterior",False)
    bm_nodes  = bm.entity_map(0).array().astype(np.int64)
    bm_coords = bm.coordinates()
    bm_cells  = bm.cells()
    bempp_boundary_grid = grid_from_element_data(bm_coords.transpose(),bm_cells.transpose())

    # Create compatible trace function space
    space = function_space(bempp_boundary_grid,"P",degree)

    # Now compute the mapping from BEM++ dofs to FEniCS dofs
    from scipy.sparse import coo_matrix
    single_order_orig = global_parameters.quadrature.far.single_order
    global_parameters.quadrature.far.single_order = 5 # TODO: use interplolation order accoriding to target space order

    def FenicsData(x, n, domain_index, result):
        value = np.empty(1)
        u_FEM.eval(value, x)
        result[0] = value

    u_FEM = dolfin.Function(fenics_space)
    N = u_FEM.vector().size()
    u_FEM.vector()[:] = np.arange(N)
    u_BEM = GridFunction(space, dual_space=space,fun=FenicsData)
    FEM2BEM = np.rint(u_BEM.coefficients).astype(np.int64)
    if max(abs(FEM2BEM-u_BEM.coefficients)) > 0.1:
        raise ValueError("interpolation from FEM to BEM space too inaccurate.")
    NB = len(FEM2BEM)
    trace_matrix = coo_matrix((
        np.ones(NB),(np.arange(NB),FEM2BEM)),
        shape=(NB,N),dtype='float64').tocsc()

    global_parameters.quadrature.far.single_order = single_order_orig

    # Now return everything
    return (space,trace_matrix)

def fenics_to_bempp_trace_data(fenics_space):
    """
    Returns tuple (space,trace_matrix)
    """

    family,degree = fenics_space_info(fenics_space)

    if family=="Lagrange":
        if degree==1:
            import p1_coupling
            return p1_coupling.p1_trace(fenics_space)
        else:
            return generic_pn_trace(fenics_space)
    else:
        raise NotImplementedError()

def fenics_space_info(fenics_space):
    """
    Returns tuple (family,degree) containing information about a FEniCS space
    """

    element = fenics_space.ufl_element()
    family = element.family()
    degree = element.degree()
    return (family,degree)
