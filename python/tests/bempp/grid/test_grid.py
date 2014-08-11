from py.test import fixture, mark


@mark.parametrize("kwargs", [
    {},
    {'topology': 'triangular'},
    {   # Incorrect topology
        'exception': ValueError,
        'topology': 'incorrect',
        'filename': "if this file exists, it's your problem"
    },
    {   # File does not exist
        'exception': IOError,
        'topology': 'triangular',
        'filename': "if this file exists, it's your problem"
    },
    {   # Cannot convert to armadillo vector of same length
        'exception': ValueError,
        'topology': 'triangular',
        'lower_left': (0., 0., 0),
        'upper_right': (1., 2.),
        'subdivisions': (4, 5)
    },
    {   # Cannot convert to armadillo vector of same type
        'exception': ValueError,
        'topology': 'triangular',
        'lower_left': ('a', 0),
        'upper_right': (1., 2.),
        'subdivisions': (4, 5)
    },
    {   # Ambiguous construction arguments:
        # Passes both filename and structured grid arguments
        'filename': "if this file exists, it's your problem",
        'topology': 'triangular',
        'lower_left': (0., 0.),
        'upper_right': (1., 2.),
        'subdivisions': (4, 5.5)
    },
    {   # Incorrect grid size
        'exception': ValueError,
        'topology': 'triangular',
        'lower_left': ('0', 0),
        'upper_right': (1., 2.),
        'subdivisions': (0, 5)
    }
])
def test_fail_on_creation(kwargs):
    """ Grid fails if arguments are incorrect """
    from py.test import raises
    from bempp.grid import Grid

    args = kwargs.copy()
    exception = args.pop('exception', TypeError)
    with raises(exception):
        return Grid(**args)


class TestGridFromMesh(object):
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
        assert grid.max_level == 0


class TestStructuredGrid(TestGridFromMesh):
    """ Creates a cartesian grid """
    @fixture
    def grid(self):
        # Fixture is also a test
        return self.test_creation()

    def test_creation(self):
        from bempp.grid import Grid
        return Grid(
            topology="triangular",
            lower_left=(0., 0.),
            upper_right=(1., 2.),
            subdivisions=(4, 5)
        )
