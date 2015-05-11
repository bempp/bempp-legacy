import numpy as _np

class FileInterface(object):

    def __init__(self, obj):
        import os.path
        from bempp import Grid
        from collections import OrderedDict

        # Reference to implementation
        self._impl = None 
        self._vertex_file_to_insertion_indices = None
        self._element_file_to_insertion_indices = None
        self._vertex_insertion_indices_to_file = []
        self._element_insertion_indices_to_file = []

        self._vertex_indices_to_file = None
        self._file_to_vertex_indices = OrderedDict()
        self._element_indices_to_file = None
        self._file_to_element_indices = OrderedDict()

        if isinstance(obj,str):
            extension = os.path.splitext(obj)[1].lower()

            if extension=='.msh':
                from bempp.file_interfaces import gmsh
                self._impl = gmsh.GmshInterface.read(obj)

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


        elif isinstance(obj,Grid):
            pass

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

    def map_vertex_to_file_key(self, vertex):
        return self._vertex_insertion_indices_to_file[self._grid.vertex_insertion_index(vertex)]

    def map_element_to_file_key(self, element):
        return self._element_insertion_indices_to_file[self._grid.element_insertion_index(element)]

    def map_element_index_to_file_key(self, index):
        return self._element_indices_to_file[index]

    def map_vertex_index_to_file_key(self, index):
        return self._vertex_indices_to_file[index]

    def map_file_key_to_element_index(self, key):
        return self._file_to_element_indices[key]

    def map_file_key_to_vertex_index(self, key):
        return self._file_to_vertex_indices[key]

    grid = property(lambda self: self._grid)


class FileInterfaceImpl(object):

    def __init__(self):
        from collections import OrderedDict

        self.__vertices = OrderedDict()
        self.__elements = OrderedDict()

    @classmethod
    def read(cls, file_name):
        pass

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





