"""Implement the mapping from FEniCS P1 to BEM++ P1 functions."""

def p1_dof_to_vertex_matrix(space):
    """Map from dofs to grid insertion vertex indices."""
    import numpy as np

    from scipy.sparse import coo_matrix

    grid = space.grid
    vertex_count = space.global_dof_count

    vertex_to_dof_map = np.zeros(vertex_count, dtype=np.int)

    for element in grid.leaf_view.entity_iterator(0):
        global_dofs = space.get_global_dofs(element)
        for ind, vertex in enumerate(element.sub_entity_iterator(2)):
            vertex_to_dof_map[grid.vertex_insertion_index(vertex)] = \
                global_dofs[ind]

    vertex_indices = np.arange(vertex_count)
    data = np.ones(vertex_count)
    return coo_matrix(
        (data, (vertex_indices, vertex_to_dof_map)), dtype='float64').tocsc()

#pylint: disable=too-many-locals
def p1_trace(fenics_space):
    """
    Return the P1 trace operator.

    This function returns a pair (space, trace_matrix),
    where space is a BEM++ space object and trace_matrix is the corresponding
    matrix that maps the coefficients of a FEniCS function to its boundary
    trace coefficients in the corresponding BEM++ space.

    """

    import dolfin #pylint: disable=import-error
    from bempp.api.fenics_interface.coupling import fenics_space_info
    from bempp.api import function_space, grid_from_element_data
    import numpy as np

    family, degree = fenics_space_info(fenics_space)
    if not (family == 'Lagrange' and degree == 1):
        raise ValueError("fenics_space must be a p1 Lagrange space")

    mesh = fenics_space.mesh()

    boundary_mesh = dolfin.BoundaryMesh(mesh, "exterior", False)
    bm_nodes = boundary_mesh.entity_map(0).array().astype(np.int64)
    bm_coords = boundary_mesh.coordinates()
    bm_cells = boundary_mesh.cells()
    bempp_boundary_grid = grid_from_element_data(
        bm_coords.transpose(), bm_cells.transpose())

    # First get trace space
    space = function_space(bempp_boundary_grid, "P", 1)

    # Now compute the mapping from FEniCS dofs to BEM++ dofs.

    # First the BEM++ dofs from the boundary vertices
    from scipy.sparse import coo_matrix
    bempp_dofs_from_b_vertices = p1_dof_to_vertex_matrix(space).transpose()

    # Now FEniCS vertices to boundary dofs
    b_vertices_from_vertices = coo_matrix(
        (np.ones(len(bm_nodes)), (np.arange(len(bm_nodes)), bm_nodes)),
        shape=(len(bm_nodes), mesh.num_vertices()),
        dtype='float64').tocsc()

    # Finally FEniCS dofs to vertices.
    vertices_from_fenics_dofs = coo_matrix(
        (np.ones(mesh.num_vertices()),
         (dolfin.dof_to_vertex_map(fenics_space), np.arange(
             mesh.num_vertices()))),
        shape=(mesh.num_vertices(), mesh.num_vertices()),
        dtype='float64').tocsc()

    # Get trace matrix by multiplication
    trace_matrix = bempp_dofs_from_b_vertices * \
        b_vertices_from_vertices * vertices_from_fenics_dofs

    # Now return everything
    return space, trace_matrix
