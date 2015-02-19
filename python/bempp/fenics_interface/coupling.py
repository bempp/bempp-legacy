import dolfin as _dolfin
from scipy.sparse import lil_matrix

def boundary_grid_from_fenics_mesh(fenics_mesh):
    """Return a boundary grid from a FEniCS Mesh."""
    try:
        fenics_mesh.bempp_boundary_grid
    except AttributeError:
        from bempp import grid_from_element_data

        bm = _dolfin.BoundaryMesh(fenics_mesh, "exterior", False)
        bm_coords = bm.coordinates()
        bm_cells  = bm.cells()
        fenics_mesh.bempp_boundary_grid = grid_from_element_data(bm_coords.transpose(),bm_cells.transpose())

    return fenics_mesh.bempp_boundary_grid

def fenics_to_bempp_map(fenics_space,bempp_space):
    """Return the permutation matrix from the FEniCS space to the BEM++ space."""
    f_to_g = fenics_to_global_map(fenics_space)
    b_to_g = bempp_to_global_map(bempp_space,fenics_space.mesh())
    return b_to_g.transpose()*f_to_g

def bempp_to_global_map(bempp_space,mesh):
    """Return the permutation matrix from the FEniCS space to the BEM++ space."""
    # This needs changing once DoF maps from BEM++ are implemented in python
    # Once changed it will not need the mesh or bd_mesh inputs.
    bd_mesh = _dolfin.BoundaryMesh(mesh, "exterior", False)
    vertices_list = bempp_space.grid.leaf_view.vertices.transpose()
    bempp_dim = bempp_space.global_dof_count
    mesh_dim = mesh.num_vertices()
    b_nodes = bd_mesh.entity_map(0).array()
    b_coords = bd_mesh.coordinates()
    not_yet_found = range(0,len(vertices_list))

    map = lil_matrix((mesh_dim,bempp_dim))
    for i,j in enumerate(b_nodes):
        for k in not_yet_found:
            if str(vertices_list[k])==str(b_coords[i]):
                map[j,k]=1
                not_yet_found.remove(k)
                break
    return map

def fenics_to_global_map(fenics_space):
    """Return the permutation matrix from the FEniCS space to the BEM++ space."""
    dim = fenics_space.dim()
    map = lil_matrix((dim,dim))
    for i,j in enumerate(_dolfin.dof_to_vertex_map(fenics_space)):
            map[j,i]=1
    return map

def fenics_to_bempp_mass_matrix(bempp_space,fenics_space):
    from bempp.operators.boundary.sparse import identity
    trace_space = boundary_space_from_fenics_space(fenics_space,bempp_space.grid)
    id = identity(trace_space,trace_space,bempp_space)
    matrix = id.weak_form().as_matrix()
    return lil_matrix(matrix)

def boundary_space_from_fenics_space(fenics_space,boundary_grid="None"):
    """Return a boundary space from a FEniCS space."""
    from bempp import function_space

    element = fenics_space.ufl_element()
    family = element.family()
    degree = element.degree()
    if boundary_grid is None:
        mesh = fenics_space.mesh()
        boundary_grid = boundary_grid_from_fenics_mesh(mesh)
    if family == "Lagrange": 
        return function_space(boundary_grid,"P",degree)
    else:
        raise NotImplementedError
