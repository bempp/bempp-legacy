#cython: embedsignature=True

<%
from data_types import dtypes, compatible_dtypes
%>

from bempp.utils cimport shared_ptr
from bempp.grid.grid cimport c_Grid, Grid
from libcpp.vector cimport vector
from bempp.utils.byte_conversion import convert_to_bytes
from bempp.utils cimport complex_float, complex_double
from cython.operator cimport dereference as deref
from bempp.assembly.grid_function cimport GridFunction
from bempp.utils.enum_types cimport gmsh_post_data_type
from cython.operator cimport dereference as deref

import os.path

__doc__="""

This module provides classes and methods to import data from Gmsh files
into BEM++ and to export BEM++ data to Gmsh.

Classes
-------

.. autoclass:: GmshInterface
    :members: write,
              save_grid_function

Functions
---------

.. autofunction:: save_grid_function_to_gmsh

"""

cdef class GmshInterface:
    """

This class provides an interface to Gmsh.
It can return a Grid object from a Gmsh file and provides conversions
between Gmsh numbering of nodes and elements and the corresponding
BEM++ numbering.

The class is initialized by either providing a valid Gmsh file
or a bempp.Grid object.

Parameters
----------
grid : bempp.Grid
    Grid to initialize the Gmsh object.
file_name : string
    Gmsh file, from which the Gmsh object is initialized.

Attributes
----------
grid: bempp.Grid
    Return a Python grid object
nodes_bempp_to_gmsh_numbering: vector[int] 
    Return a vector v, where v[i] gives the Gmsh number of the ith node in BEM++.
elements_bempp_to_gmsh_numbering: vector[int]
    Return a vector v, where v[i] gives the Gmsh number of the ith element in BEM++.
nodes_gmsh_to_bempp_numbering: vector[int]
    Return a vector v, where v[i] gives the BEM++ number of the node with number i in Gmsh 
    or -1 if the node is not contained in the BEM++ grid object.
elements_gmsh_to_bempp_numbering: vector[int]
    Return a vector v, where v[i] gives the BEM++ number of the node with number i in Gmsh
    or -1 if the node is not contained in the BEM++ grid object.


Examples
--------
To load the grid from a Gmsh file named ``grid.msh`` use

>>> grid = GmshInterface(file_name="grid.msh").grid()

To write an existing grid object to a new file ``grid.msh`` use

>>> GmshInterface(grid=grid).write("grid.msh")

    """


    def __cinit__(self,grid=None,file_name=None):
        pass

    def __init__(self,grid=None,file_name=None): 


        if ((grid is not None and file_name is not None) or
                (grid is None and file_name is None)):
            raise ValueError("Exactly one of `file_name` or `grid` must be specified")

        physical_entity = -1 # At the moment read all physical entities by default

        if file_name is not None:

            if not os.path.isfile(file_name):
                raise ValueError("File does not exist")
            self.impl_.reset(new c_GmshIo(convert_to_bytes(file_name),physical_entity))

        if grid is not None:

            if not isinstance(grid, Grid):
                raise ValueError("grid is not a valid Grid object")

            self.impl_.reset(new c_GmshIo((<Grid>grid).impl_))

    property grid:
        """ Return a Grid object. """

        def __get__(self):

            cdef Grid grid = Grid.__new__(Grid)
            grid.impl_ = deref(self.impl_).grid()
            return grid

    property nodes_bempp_to_gmsh_numbering:
        """ Return a vector v, where v[i] gives the Gmsh number of ith node in BEM++. """

        def __get__(self):

            return deref(self.impl_).nodePermutation()


    property elements_bempp_to_gmsh_numbering:
        """ Return a vector v, where v[i] gives the Gmsh number of ith element in BEM++. """

        def __get__(self):

            return deref(self.impl_).elementPermutation()
        

    property nodes_gmsh_to_bempp_numbering:
        """ Return a vector v, where v[i] gives the BEM++ number of the node with number i in Gmsh. 

            Notes
            -----
            If v[i] = -1 the corresponding element is not contained in the BEM++ grid object.

        """

        def __get__(self):

            return deref(self.impl_).inverseNodePermutation()

    property elements_gmsh_to_bempp_numbering:
        """ Return a vector v, where v[i] gives the BEM++ number of the element with number i in Gmsh. 
            
            Notes
            -----
            If v[i] = -1 the corresponding node is not contained in the BEM++ grid object.
        
        """

        def __get__(self):

            return deref(self.impl_).inverseElementPermutation()

    def write(self,file_name):
        """ 
        
        Write data in Gmsh format. 

        Parameters
        ----------
        file_name : string
            Name of file to write to (file is overwritten).

        """

        deref(self.impl_).write(convert_to_bytes(file_name))

    def save_grid_function(self,GridFunction grid_function, object data_label, object gmsh_type="element_node",
            object complex_mode="real"):
        """ 
        
        Store a GridFunction in the Gmsh object


        Parameters
        ----------
        grid_function : bempp.GridFunction
            Gridfunction object to save.
        data_label : string
            Name of data set.
        gmsh_type : string
            Type of gmsh post processing data to save. Options are
            "element", "node" or "element_node".
        complex_mode : string
            Determines what part of a complex data set is stored. possible
            options are "real", "imag", "abs" or "all". If "all" is chosen
            then different data sets are saved, in the same file with
            the real part, imaginary part, and absolute value of the data.


        See Also
        --------
        save_grid_function_to_gmsh

        Examples
        --------
        To store a GridFunction object grid_function use the command

        >>> gmsh_interface.save_grid_function(grid_function,"data")
        >>> gmsh_interface.write("data.msh")

        """

% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]: 

        if grid_function.basis_type == "${pybasis}" and grid_function.result_type=="${pyresult}":
            c_exportToGmsh[${cybasis},${cyresult}](
                    deref(grid_function._impl_${pybasis}_${pyresult}),
                    convert_to_bytes(data_label),deref(self.impl_), 
                    gmsh_post_data_type(convert_to_bytes(gmsh_type)),convert_to_bytes(complex_mode))
            return
%         endif
%     endfor
% endfor


def save_grid_function_to_gmsh(grid_function, data_label, file_name, 
        gmsh_type="element_node", complex_mode="real"):
   """ 
   
   Save a GridFunction to a Gmsh file.

   Parameters
   ----------
   grid_function : bempp.GridFunction
       Gridfunction object to save.
   data_label : string
       Name of data set.
   file_name : string
       Name of output file.
   gmsh_type : string
        Type of gmsh post processing data to save. Options are
        "element", "node" or "element_node".
   complex_mode : string
       Determines what part of a complex data set is stored. possible
       options are "real", "imag", "abs" or "all". If "all" is chosen
       then different data sets are saved, in the same file with
       the real part, imaginary part, and absolute value of the data.

   See Also
   --------
   GmshInterface.save_grid_function

   Examples
   --------
   To store a GridFunction object grid_function use the command

   >>> save_grid_function_to_gmsh(grid_function,"data","data.msh")

   """
   gmsh_interface = GmshInterface(grid=grid_function.grid)
   gmsh_interface.save_grid_function(<GridFunction>grid_function,data_label,gmsh_type,complex_mode)
   gmsh_interface.write(file_name)





