import numpy as _np

class FileReader(object):
    """

    This class represents a generic interface to read data from different
    file sources. Currently it supports the following types:

    * Gmsh ASCII v2.2 files

    The class creates a grid object from the data and provides maps
    that translate between the internal BEM++ numbering for elements
    and vertices and the file numbering.

    All parameters for the constructor are keyword arguments.

    Parameters
    ----------
    file_name : string
        Name of the input file.

    Attributes
    ----------
    grid : bempp.api.Grid
        Returns a grid object, representing the element and node data
        in the file.
    vertex_index_to_file_key_map : list
        vertex_index_to_file_key_map[i] returns the associated file key
        of the ith vertex in BEM++.
    vertex_file_key_to_index_map : dictionary
        vertex_file_key_to_index_map[key] returns the BEM++ index of the
        vertex with name 'key' in the file.
    element_index_to_file_key_map : list
        element_index_to_file_key_map[i] returns the associated file key
        of the ith element in BEM++.
    element_file_key_to_index_map : dictionary
        element_file_key_to_index_map[key] returns the BEM++ index of the
        element with name 'key' in the file.

    """


    def __init__(self,**kwargs):
        import os.path
        from collections import OrderedDict

        self._impl = None 
        self._vertex_file_to_insertion_indices = None
        self._element_file_to_insertion_indices = None
        self._vertex_insertion_indices_to_file = []
        self._element_insertion_indices_to_file = []

        self._vertex_indices_to_file = None
        self._file_to_vertex_indices = OrderedDict()
        self._element_indices_to_file = None
        self._file_to_element_indices = OrderedDict()

        if 'file_name' in kwargs:
            fname = kwargs['file_name']
            extension = os.path.splitext(fname)[1].lower()

            if extension=='.msh':
                from bempp.api.file_interfaces import gmsh
                self._impl = gmsh.GmshInterface.read(fname)

            self._grid = self._create_grid()

            # Setup the mappers

            self._file_to_vertex_indices = OrderedDict.fromkeys(self._impl.vertices.keys(),value=-1)
            self._file_to_element_indices = OrderedDict.fromkeys(self._impl.elements.keys(),value=-1)

            self._vertex_indices_to_file = self._grid.leaf_view.entity_count(2)*[None]
            self._element_indices_to_file = self._grid.leaf_view.entity_count(0)*[None]

            index_set = self._grid.leaf_view.index_set()
            for elem in self._grid.leaf_view.entity_iterator(0):
                index = index_set.entity_index(elem)
                insertion_index = self._grid.element_insertion_index(elem)
                file_key = self._element_insertion_indices_to_file[insertion_index]
                self._element_indices_to_file[index] = file_key
                self._file_to_element_indices[file_key] = index
                
            for vertex in self._grid.leaf_view.entity_iterator(2):
                index = index_set.entity_index(vertex)
                insertion_index = self._grid.vertex_insertion_index(vertex)
                file_key = self._vertex_insertion_indices_to_file[insertion_index]
                self._vertex_indices_to_file[index] = file_key
                self._file_to_vertex_indices[file_key] = index

    def _create_grid(self):

        from bempp.api import GridFactory
        from collections import OrderedDict

        vertices = self._impl.vertices
        elements = self._impl.elements

        factory = GridFactory()

        self._element_file_to_insertion_indices = OrderedDict.fromkeys(elements.keys(),value=-1)
        self._vertex_file_to_insertion_indices = OrderedDict.fromkeys(vertices.keys(),value=-1)

        vertex_insertion_count = 0
        element_insertion_count = 0

        for key in elements:
            elem = elements[key]
            elem_vertex_keys = []
            for vertex_key in elem['data']:
                if self._vertex_file_to_insertion_indices[vertex_key] ==-1:
                    factory.insert_vertex(vertices[vertex_key])
                    self._vertex_file_to_insertion_indices[vertex_key] = vertex_insertion_count
                    self._vertex_insertion_indices_to_file.append(vertex_key)
                    vertex_insertion_count += 1
                elem_vertex_keys.append(self._vertex_file_to_insertion_indices[vertex_key])
            factory.insert_element(elem_vertex_keys, domain_index=elem['domain_index'])
            self._element_file_to_insertion_indices[key] = element_insertion_count
            self._element_insertion_indices_to_file.append(key)
            element_insertion_count += 1

        return factory.finalize()

    vertex_index_to_file_key_map = property(lambda self: self._vertex_indices_to_file)
    vertex_file_key_to_index_map = property(lambda self: self._file_to_vertex_indices)

    element_index_to_file_key_map = property(lambda self: self._element_indices_to_file)
    element_file_key_to_index_map = property(lambda self: self._file_to_element_indices)

    grid = property(lambda self: self._grid)

def import_grid(file_name):
    """

    A simple grid importer. Use this instead of FileReader if
    information about the numbering of vertices and elements
    from the grid file needs not be preserved.

    Parameters
    ----------
    file_name : string
       Name of the file from which to read the grid.

    Returns
    -------
    grid : bempp.api.Grid
        A grid object

    Examples
    --------
    To read a Gmsh grid file called 'grid_file.msh' use

    >>> grid = import_grid('grid_file.msh')

    """

    return FileReader(file_name=file_name).grid


def export(**kwargs):
    """

    This funtion can export grids and gridfunctions into external file formats.
    Supported formats are:

    * Gmsh ASCII v2.2 files

    The function only takes keyword arguments. 

    Parameters
    ----------
    file_name : string
        Name of the output file.
    grid : bempp.api.Grid
        A grid object to write out
    grid_function : bempp.api.GridFunction
        A gridfunction to write out
    vertex_index_to_file_key_map :  list
        (optional) A list that maps BEM++ indices to vertex indices in the file
    element_index_to_file_key_map : list
        (optional) A list that maps BEM++ indices to element indices in the file
    data_type : string
        (optional) One of 'node', 'element' or 'element_node' and specifies if
        data in grid functions is associated with nodes, elements or elementwise
        nodes in the data file. The default behavior is given by the specific
        file writer for the datatype (e.g. 'element_node' for Gmsh) but can
        be overridden here.
    label : string
        (optional) A string labelling grid function data in the file.
    transformation : function object
        (optional) A function object that is applied to the data before
        writing it out.

    Notes
    -----
    * A grid or grid_function object must always be specified.
    * Depending on the file format an offset of 1 is added to the indices if
      no index_to_file_key_map is provided. The reason is that some
      file formats such as Gmsh start counting nodes and elements from 1
      instead of zero.

    Examples
    --------
    To save the grid object 'mygrid' in Gmsh format use

    >>> export(grid=mygrid, file_name='output.msh')

    To save the grid_function object 'gridfun' in Gmsh format
    with data associated to elements use

    >>> export(grid_function=gridfun, file_name='output.msh', data_type='element')

    To save the real part of a complex grid function in Gmsh format use

    >>> export(grid_function=gridfun, file_name='output.msh', transformation=lambda x: np.real(x))

    """

    import os

    interface = None # Holds the actual FileInterface for the specified data format
    vertex_index_to_file_key_map = None
    element_index_to_file_key_map = None

    if 'file_name' in kwargs:
        fname = kwargs['file_name']
    else:
        raise ValueError("file_name must be specified.")
    
    extension = os.path.splitext(fname)[1].lower()

    if extension=='.msh':
        from bempp.api.file_interfaces import gmsh
        interface = gmsh.GmshInterface()
    
    if int('grid' in kwargs) + int('grid_function' in kwargs) != 1:
        raise ValueError("Exactly one of 'grid' or 'grid_function' must be specified")

    if 'grid' in kwargs:
        grid = kwargs['grid']
    elif 'grid_function' in kwargs:
        grid = kwargs['grid_function'].grid

    number_of_vertices = grid.leaf_view.entity_count(2)
    number_of_elements = grid.leaf_view.entity_count(0)

    offset = interface.index_offset
    if 'vertex_index_to_file_key_map' in kwargs:
        vertex_index_to_file_key_map = kwargs['vertex_index_to_file_key_map']
    else:
        vertex_index_to_file_key_map = range(offset,number_of_vertices+offset)
    if 'element_index_to_file_key_map' in kwargs:
        element_index_to_file_key_map = kwargs['element_index_to_file_key_map']
    else:
        element_index_to_file_key_map = range(offset,number_of_elements+offset)

    # Create the vertex and element structure

    from collections import OrderedDict

    vertex_iterator = grid.leaf_view.entity_iterator(2)
    element_iterator = grid.leaf_view.entity_iterator(0)
    index_set = grid.leaf_view.index_set()

    vertices = OrderedDict([(vertex_index_to_file_key_map[index_set.entity_index(vertex)],vertex.geometry.corners[:,0])
        for vertex in vertex_iterator])
    elements = OrderedDict([(element_index_to_file_key_map[index_set.entity_index(element)],
        {'data':[vertex_index_to_file_key_map[index_set.sub_entity_index(element,n,2)] for n in range(3)],
         'domain_index':element.domain}) for element in element_iterator])

    interface.add_grid_data(vertices,elements)

    # Evaluate data

    if 'grid_function' in kwargs:
        fun = kwargs['grid_function']
        data_type = kwargs.get('data_type',interface.default_data_type)

        if 'transformation' in kwargs:
            transformation = kwargs['transformation']
        else:
            transformation = lambda x: x

        index_set = grid.leaf_view.index_set()

        if data_type == 'element_node':
            local_coordinates = _np.array([[0,1,0],[0,0,1]])
            data = OrderedDict.fromkeys(element_index_to_file_key_map)

            for element in grid.leaf_view.entity_iterator(0):
                data[element_index_to_file_key_map[index_set.entity_index(element)]] = transformation(
                        fun.evaluate(element,local_coordinates))
            interface.add_element_node_data(data,kwargs.get('label','element_node_data'))
        elif data_type == 'node':
            local_coordinates = _np.array([[0,1,0],[0,0,1]])
            data = OrderedDict.fromkeys(vertex_index_to_file_key_map)
            for element in grid.leaf_view.entity_iterator(0):
                local_data = transformation(fun.evaluate(element,local_coordinates))
                for i in range(3):
                    data[vertex_index_to_file_key_map[index_set.sub_entity_index(element,i,2)]] = local_data[:,i]
            interface.add_node_data(data,kwargs.get('label','node_data'))
        elif data_type == 'element':
            local_coordinates = _np.array([[1./3],[1./3]])
            data = OrderedDict.fromkeys(element_index_to_file_key_map)

            for element in grid.leaf_view.entity_iterator(0):
                data[element_index_to_file_key_map[index_set.entity_index(element)]] = transformation(
                        fun.evaluate(element,local_coordinates).ravel())
            interface.add_element_data(data,kwargs.get('label','element_data'))
        else:
            raise ValueError("data_type must be one of 'node', 'element', or 'element_node'")

    interface.write(kwargs['file_name'])


class FileInterfaceImpl(object):

    def __init__(self):
        from collections import OrderedDict

        self.__vertices = OrderedDict()
        self.__elements = OrderedDict()

    @classmethod
    def read(cls, file_name):
        pass

    def write(self, file_name):
        pass

    def add_grid_data(self, vertices, elements):
        self.__vertices = vertices
        self.__elements = elements

    def add_node_data(self, data, label):
        pass

    def add_element_data(self, data, label):
        pass

    def add_element_node_data(self, data, label):
        pass

    @property
    def index_offset(self):
        return 1

    @property
    def default_data_type(self):
        return 'node'
        

    vertices = property(lambda self: self.__vertices)
    elements = property(lambda self: self.__elements)
    

class Vertex(object):

    def __init__(self, index, x, y, z):
        self.index = index
        self.data = _np.array([x,y,z],dtype='float64')

class Element(object):

    def __init__(self, index, v0, v1, v2, domain_index = 0):
        self.index = index
        self.data = [v0, v1, v2]
        self.domain_index = domain_index
