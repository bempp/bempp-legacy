from py.test import fixture, mark


class TestStructuredGrid(object):
    """ Creates a cartesian grid """

    def test_creation(self):
        from bempp.grid import structured_grid
        return structured_grid(
            (0., 0.),
            ('1.', '2.'),
            (4, 5)
        )


class TestConnectivityGrid(object):
    """ Creates grid from connectivity data """

    vertices = [[0, '0', 1, 1.2], [0, 1, 0, 1.1], [0, 0, 0, 0.5]]
    corners = [[0, 2], [1, 1], [2, 3], [-1, -1]]

    def test_creation(self):
        from bempp.grid import grid_from_element_data
        return grid_from_element_data(
                self.vertices,self.corners
        )

