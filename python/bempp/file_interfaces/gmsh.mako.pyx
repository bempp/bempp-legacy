#cython: embedsignature=True


from bempp.utils cimport shared_ptr
from bempp.grid.grid cimport c_Grid, Grid
from libcpp.vector cimport vector
from bempp.utils.byte_conversion import convert_to_bytes
from cython.operator cimport dereference as deref

import os.path

__doc__="""

This module provides classes and methods to import data from Gmsh files
into BEM++ and to export BEM++ data to Gmsh.

Classes
-------

.. autoclass:: Gmsh
    :members: write

"""

cdef class Gmsh:
    """

This class provides an interface to an existing Gmsh .msh file.
It can return a Grid object from a Gmsh file and provides conversions
between Gmsh numbering of nodes and elements and the corresponding
BEM++ numbering.

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


Example
-------
To load the grid from a Gmsh file named ``grid.msh`` use

>>> grid = Gmsh("grid.msh").grid()

    """


    def __cinit__(self,**kwargs):
        pass

    def __init__(self,file_name): 

        physical_entity = -1 # At the moment read all physical entities by default

        if not os.path.isfile(file_name):
            raise ValueError("File does not exist")
            self.impl_.reset(new c_GmshIo(convert_to_bytes(kwargs['file_name']),physical_entity))

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
        """ Write out grid data in Gmsh format. 

            Parameters
            ----------
            file_name : string
                name of file to write to (file is overwritten)

        """

        deref(self.impl_).write(convert_to_bytes(file_name))

