from bempp.utils cimport shared_ptr
from bempp.grid.grid cimport c_Grid, Grid
from libcpp.vector cimport vector
from bempp.utils.byte_conversion import convert_to_bytes
from cython.operator cimport dereference as deref

import os.path

cdef class Gmsh:

    def __cinit__(self, grid=None, file_name=None, physical_entity=-1):
        pass

    def __init__(self, grid=None, file_name=None, physical_entity=-1):

        if grid is not None:
            self.impl_.reset(new c_GmshIo((<Grid>grid).impl_))
        elif file_name is not None:
            if not os.path.isfile(file_name):
                raise ValueError("File does not exist")
            self.impl_.reset(new c_GmshIo(convert_to_bytes(file_name),physical_entity))
        else:
            raise ValueError("One of `grid` or `file_name` must be provided.")


    property grid:
        """ Return a bempp.grid object. """

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
        """ gmsh.write(file_name) 
        
            Write out a grid in Gmsh format with name given by 'file_name'. 

            Parameters
            ----------
            
            file_name : string
                Name of output file

        """

        deref(self.impl_).write(convert_to_bytes(file_name))


