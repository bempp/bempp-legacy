def p1_dof_to_vertex_matrix(space):
    """Return permutation matrix that maps from dofs to grid insertion vertex indices."""
    import numpy as np

    from scipy.sparse import coo_matrix

    grid = space.grid
    vertex_count = space.global_dof_count

    vertex_to_dof_map = np.zeros(vertex_count, dtype=np.int)

    for element in grid.leaf_view.entity_iterator(0):
        global_dofs = space.get_global_dofs(element)
        for ind, v in enumerate(element.sub_entity_iterator(2)):
            vertex_to_dof_map[grid.vertex_insertion_index(v)] = global_dofs[ind]

    vertex_indices = np.arange(vertex_count)
    data = np.ones(vertex_count)
    return coo_matrix((data, (vertex_indices, vertex_to_dof_map)), dtype='float64').tocsc()


def p1_trace(fenics_space):
    import dolfin
    from .coupling import fenics_space_info
    from bempp import function_space, grid_from_element_data
    import numpy as np

    family, degree = fenics_space_info(fenics_space)
    if not (family == 'Lagrange' and degree == 1):
        raise ValueError("fenics_space must be a p1 Lagrange space")

    mesh = fenics_space.mesh()

    bm = dolfin.BoundaryMesh(mesh, "exterior", False)
    bm_nodes = bm.entity_map(0).array().astype(np.int64)
    bm_coords = bm.coordinates()
    bm_cells = bm.cells()
    bempp_boundary_grid = grid_from_element_data(bm_coords.transpose(), bm_cells.transpose())

    # First get trace space 
    space = function_space(bempp_boundary_grid, "P", 1)

    # Now compute the mapping from BEM++ dofs to FEniCS dofs

    # First the BEM++ dofs to the boundary vertices
    from scipy.sparse import coo_matrix

    vertex_indices = np.arange(space.global_dof_count)
    data = np.ones(space.global_dof_count)
    bempp_dofs_from_b_vertices = p1_dof_to_vertex_matrix(space).transpose()

    # Now the boundary vertices to all the vertices
    b_vertices_from_vertices = coo_matrix((
        np.ones(len(bm_nodes)), (np.arange(len(bm_nodes)), bm_nodes)),
        shape=(len(bm_nodes), mesh.num_vertices()), dtype='float64').tocsc()

    # Finally the vertices to FEniCS dofs
    vertices_from_fenics_dofs = coo_matrix((
        np.ones(mesh.num_vertices()), (dolfin.dof_to_vertex_map(fenics_space), np.arange(mesh.num_vertices()))),
        shape=(mesh.num_vertices(), mesh.num_vertices()), dtype='float64').tocsc()

    # Get trace matrix by multiplication
    trace_matrix = bempp_dofs_from_b_vertices * b_vertices_from_vertices * vertices_from_fenics_dofs

    # Now return everything
    return space, trace_matrix
