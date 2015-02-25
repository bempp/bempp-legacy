
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

    # First the BEM++ dofs to the boundary vertices
    from ._lagrange_coupling import p1_vertex_map
    from scipy.sparse import coo_matrix
    # this is giving [0,1,2,...]
    vertex_to_dof_map =  p1_vertex_map(space)
    # this is giving [0,1,2,...]
    vertex_indices = np.arange(space.global_dof_count)
    data = np.ones(space.global_dof_count)
    # this gives identity matrix...
    bempp_dofs_from_b_vertices = coo_matrix((data,(vertex_to_dof_map,vertex_indices)),dtype='float64').tocsr()
    # this gives identity matrix...

    # Now the boundary vertices to all the vertices
    b_vertices_from_vertices = coo_matrix((
        np.ones(len(bm_nodes)),(np.arange(len(bm_nodes)),bm_nodes)),
        shape=(len(bm_nodes),mesh.num_vertices()),dtype='float64').tocsr()

    # Finally the vertices to FEniCS dofs
    vertices_from_fenics_dofs = coo_matrix((
        np.ones(mesh.num_vertices()),(dolfin.dof_to_vertex_map(fenics_space),np.arange(mesh.num_vertices()))),
        shape=(mesh.num_vertices(),mesh.num_vertices()),dtype='float64').tocsr()

    # Get trace matrix by multiplication
    trace_matrix = bempp_dofs_from_b_vertices*b_vertices_from_vertices*vertices_from_fenics_dofs

    # Now compute the mass matrix 
    from bempp.operators.boundary.sparse import identity
    mass_matrix = trace_matrix.transpose()*identity(space,space,space).weak_form().sparse_operator

    # Now return everything
    return (space,mass_matrix,trace_matrix,bempp_dofs_from_b_vertices)
