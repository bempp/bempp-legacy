from py.test import mark, fixture
from numpy import dtype
from bempp.file_interfaces import gmsh

_mesh_name = "sphere-h-0.4.msh"

class TestGmsh(object):


    @fixture
    def mesh_path(self):
        from os.path import join, exists
        from bempp.config import paths
        filename = join(paths.meshes, _mesh_name)
        if not exists(filename):
            raise IOError("Mesh %s does not exist" % filename)
        return filename

    @fixture
    def grid(self):
        # Fixture is also a test
        return self.test_instantiate_from_file()

    def test_instantiate_from_grid(self):
        gmsh_interface = gmsh.GmshInterface(self.grid())

    def test_instantiate_from_file(self):
        gmsh_interface = gmsh.GmshInterface(self.mesh_path())
        return gmsh_interface.grid
