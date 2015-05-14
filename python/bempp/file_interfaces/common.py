import numpy as _np

class FileReader(object):

    def __init__(self,**kwargs):
        import os.path
        from bempp import Grid
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

        if kwargs.has_key('file_name'):
            fname = kwargs['file_name']
            extension = os.path.splitext(fname)[1].lower()

            if extension=='.msh':
                from bempp.file_interfaces import gmsh
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

        from bempp import GridFactory
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


def export(**kwargs):

        import os

        interface = None # Holds the actual FileInterface for the specified data format
        vertex_index_to_file_key_map = None
        element_index_to_file_key_map = None

        if kwargs.has_key('file_name'):
            fname = kwargs['file_name']
        else:
            raise ValueError("file_name must be specified.")
        
        extension = os.path.splitext(fname)[1].lower()

        if extension=='.msh':
            from bempp.file_interfaces import gmsh
            interface = gmsh.GmshInterface()
        
        if kwargs.has_key('grid') + kwargs.has_key('grid_function')!= 1:
            raise ValueError("Exactly one of 'grid' or 'grid_function' must be specified")

        if kwargs.has_key('grid'):
            grid = kwargs['grid']
        elif kwargs.has_key('grid_function'):
            grid = kwargs['grid_function'].grid

        number_of_vertices = grid.leaf_view.entity_count(2)
        number_of_elements = grid.leaf_view.entity_count(0)

        if kwargs.has_key('vertex_index_to_file_key_map'):
            vertex_index_to_file_key_map = kwargs['vertex_index_to_file_key_map']
        else:
            vertex_index_to_file_key_map = range(number_of_vertices)
        if kwargs.has_key('element_index_to_file_key_map'):
            element_index_to_file_key_map = kwargs['element_index_to_file_key_map']
        else:
            element_index_to_file_key_map = range(number_of_elements)

        # Create the vertex and element structure

        from collections import OrderedDict

        vertex_iterator = grid.leaf_view.entity_iterator(2)
        element_iterator = grid.leaf_view.entity_iterator(0)
        index_set = grid.leaf_view.index_set()

        vertices = OrderedDict([(vertex_index_to_file_key_map[index_set.entity_index(vertex)],vertex.geometry.corners[:,0])
            for vertex in vertex_iterator])
        elements = OrderedDict([(element_index_to_file_key_map[index_set.entity_index(element)],
            {'data':[element_index_to_file_key_map[index_set.sub_entity_index(element,n,2)] for n in range(3)],
             'domain_index':element.domain}) for element in element_iterator])

        interface.add_grid_data(vertices,elements)

        # Evaluate data

        if kwargs.has_key('grid_function'):
            fun = kwargs['grid_function']
            data_type = kwargs.get('data_type',interface.default_data_type)

            if kwargs.has_key('transformation'):
                transformation = kwargs['transformation']
            else:
                transformation = lambda x: x

            index_set = grid.leaf_view.index_set()

            if data_type == 'element_node':
                local_coordinates = _np.array([[0,1,0],[0,0,1]])
                data = OrderedDict.fromkeys(element_index_to_file_key_map)

                for element in grid.leaf_view.entity_iterator(0):
                    data[element_index_to_file_key_map[index_set.entity_index(element)]] = transformation(
                            fun.evaluate(element,local_coordinates)).ravel()
                interface.add_data_set(data,'element_node',kwargs.get('label','element_node_data'))
            elif data_type == 'node':
                local_coordinates = _np.array([[0,1,0],[0,0,1]])
                data = OrderedDict.fromkeys(vertex_index_to_file_key_map)
                for element in grid.leaf_view.entity_iterator(0):
                    local_data = transformation(fun.evaluate(element,local_coordinates)).ravel()
                    for i in range(3):
                        data[vertex_index_to_file_key_map[index_set.sub_entity_index(element,i,2)]] = local_data[i]
                interface.add_data_set(data,'node',kwargs.get('label','node_data'))
            elif data_type == 'element':
                local_coordinates = _np.array([[1./3],[1./3]])
                data = OrderedDict.fromkeys(element_index_to_file_key_map)

                for element in grid.leaf_view.entity_iterator(0):
                    data[element_index_to_file_key_map[index_set.entity_index(element)]] = transformation(
                            fun.evaluate(element,local_coordinates)).ravel()
                interface.add_data_set(data,'element',kwargs.get('label','element_data'))
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

    def add_data_set(self, data, data_type, label):
        pass

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





