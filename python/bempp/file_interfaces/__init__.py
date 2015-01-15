from bempp.grid import Grid
from bempp import GridFunction


def export(obj, file_name,**kwargs):
    import os.path

    extension = os.path.splitext(file_name)[1].lower()

    if extension=='.msh':
        from bempp.file_interfaces import gmsh

        if isinstance(obj,Grid):
            gmsh.GmshInterface(obj).write(file_name)
        elif isinstance(obj,GridFunction):
            data_label = kwargs['data_label']
            del kwargs['data_label']
            gmsh.save_grid_function_to_gmsh(obj, data_label, file_name, **kwargs)

    else:
        raise ValueError("Unknown file extension")

        
        




