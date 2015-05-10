import numpy as _np

class FileInterface(object):

    def __init__(self):
        self._nodes = None
        self._elements = None

    def _create_grid(self):


    @classmethod
    def read(cls, file_name):
        pass

    def vertex_bempp_to_file_index(self, bempp_index):
        pass

    def vertex_file_to_bempp_index(self, file_index):
        pass

    def element_bempp_to_file_index(self, bempp_index):
        pass

    def element_file_to_bempp_index(self, file_index):
        pass


class Vertex(object):

    def __init__(self, index, x, y, z):
        self.index = index
        self.data = _np.array([x,y,z],dtype='float64')

class Element(object):

    def __init__(self, index, v0, v1, v2, domain_index = 0):
        self.index = index
        self.data = [v0, v1, v2]
        self.domain_index = domain_index





