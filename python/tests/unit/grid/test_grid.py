from py.test import fixture


class TestGridFromMesh:
    """ Load grid from file """
    @fixture
    def mesh_path(self):
        from os.path import join, exists
        from bempp.config import paths
        filename = join(paths.meshes, "sphere-h-0.4.msh")
        if not exists(filename):
            raise IOError("Mesh %s does not exist" % filename)
        return filename

    @fixture
    def grid(self, mesh_path):
        # Fixture is also a test
        return self.test_creation(mesh_path)

    def test_creation(self, mesh_path):
        from bempp.grid import Grid
        return Grid(topology="triangular", filename=mesh_path)

    def test_properties(self, grid):
        assert grid.topology == "triangular"
        assert grid.dim == 2
        assert grid.dim_world == 3
