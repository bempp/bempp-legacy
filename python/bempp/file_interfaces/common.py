import numpy as _np

class FileInterface(object):

    def write(self, file_name):
        pass

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

def node(index, x, y, z):

    return {'index': index, 'data':_np.array([x,y,z],dtype='float64')}

def element(index, v0, v1, v2):

    return {'index': index,
            'data': [v0, v1, v2]}






