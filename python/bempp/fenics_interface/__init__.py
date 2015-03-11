try:
    import dolfin as _dolfin
except:
    print("Cannot find dolfin. Please check your FEniCS installation.")

from .fenics_operator import FenicsOperator
from .coupling import boundary_grid_from_fenics_mesh,fenics_to_bempp_trace_data,fenics_space_info
