import dolfin as _dolfin


def boundary_grid_from_fenics_mesh(fenics_mesh):
    """

    Return a boundary grid from a FEniCS Mesh

    """
    from bempp import grid_from_element_data
    

    bm = _dolfin.BoundaryMesh(fenics_mesh, "exterior", False)
    bm_coords = bm.coordinates()
    bm_cells  = bm.cells()
    grid = grid_from_element_data(bm_coords.transpose(),bm_cells.transpose())

    return grid

def fenics_to_bempp_coord_transformation(bempp_space,fenics_space):
    pass

def fenics_to_bempp_mass_matrix(bempp_space,fenics_space):
    pass

def boundary_space_from_fenics_space(fenics_space):
    pass

def _boundary_space_from_p1(fenics_space):
    
    mesh = fenics_space.mesh()

