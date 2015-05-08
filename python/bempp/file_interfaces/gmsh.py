from .common import FileInterface


class GmshInterface(FileInterface):

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

