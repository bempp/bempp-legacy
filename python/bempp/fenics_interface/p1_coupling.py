
def p1_dof_to_vertex_matrix(space):

    import numpy as np

    from ._lagrange_coupling import p1_vertex_map
    from scipy.sparse import coo_matrix
    vertex_to_dof_map =  p1_vertex_map(space)
    vertex_indices = np.arange(space.global_dof_count)
    data = np.ones(space.global_dof_count)
    return coo_matrix((data,(vertex_indices,vertex_to_dof_map)),dtype='float64').tocsr()
    
def p1_trace(fenics_space):

    import dolfin
    from .coupling import fenics_space_info
    from bempp import function_space,grid_from_element_data
    import numpy as np

    family,degree = fenics_space_info(fenics_space)
    if not (family=='Lagrange' and degree == 1):
        raise ValueError("fenics_space must be a p1 Lagrange space")

    mesh = fenics_space.mesh()

    bm = dolfin.BoundaryMesh(mesh,"exterior",False)
    bm_nodes  = bm.entity_map(0).array().astype(np.int64)
    bm_coords = bm.coordinates()
    bm_cells  = bm.cells()
    bempp_boundary_grid = grid_from_element_data(bm_coords.transpose(),bm_cells.transpose())

    # First get trace space 

    space = function_space(bempp_boundary_grid,"P",1)

    # Now compute the mapping from BEM++ dofs to FEniCS dofs

    # First the BEM++ dofs to the vertices

    from ._lagrange_coupling import p1_vertex_map
    from scipy.sparse import coo_matrix
    vertex_to_dof_map =  p1_vertex_map(space)
    vertex_indices = np.arange(space.global_dof_count)
    data = np.ones(space.global_dof_count)
    bempp_dofs_to_vertices = coo_matrix((data,(vertex_indices,vertex_to_dof_map)),dtype='float64').tocsr()

    # Now the vertices to FEniCS dofs

    vertices_to_fenics_dofs = coo_matrix((
        np.ones(len(bm_nodes)),(bm_nodes,np.arange(len(bm_nodes)))),
        shape=(mesh.num_vertices(),len(bm_nodes)),dtype='float64').tocsr()

    # Get trace matrix by multiplication
    # FIXME Remove transpose

    trace_matrix = bempp_dofs_to_vertices.transpose()*vertices_to_fenics_dofs.transpose()
    
    # Now compute the mass matrix 

    from bempp.operators.boundary.sparse import identity

    mass_matrix = trace_matrix.transpose()*identity(space,space,space).weak_form().sparse_operator

    # Now return everything

    return (space,mass_matrix,trace_matrix)











