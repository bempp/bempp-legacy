"""Unit tests for the grid class."""
import unittest
import numpy as np

VERTICES = np.array([[0, 1, 1, 0],
                     [0, 0, 1, 1],
                     [0, 0, 0, 0]])

ELEMENTS = np.array([[0, 1],
                     [1, 2],
                     [3, 3]])


class TestGrid(unittest.TestCase):
    """Test the Grid implementation."""

    def test_structured_grid(self):
        """Test structured grid generation."""
        from bempp.api import structured_grid

        n_x = 4
        n_y = 5

        grid = structured_grid(
            (0., 0.),
            ('1.', '2.'),
            (n_x, n_y)
        )

        box = grid.bounding_box

        self.assertTrue(np.all(box == np.array([[0, 0, 0], [1, 2, 0]])))
        self.assertEqual(grid.leaf_view.entity_count(0), 2 * n_x * n_y)

    def test_grid_from_element_data(self):
        """Test grid generation from element data."""
        from bempp.api.grid.grid import grid_from_element_data

        grid = grid_from_element_data(VERTICES, ELEMENTS)

        box = grid.bounding_box
        self.assertTrue(np.all(box == np.array([[0, 0, 0], [1, 1, 0]])))
        self.assertEqual(grid.leaf_view.entity_count(0), 2)
        self.assertEqual(grid.leaf_view.entity_count(2), 4)

        actual_vertices = grid.leaf_view.vertices
        actual_elements = grid.leaf_view.elements

        self.assertEqual(actual_vertices.shape[1], 4)
        self.assertEqual(actual_elements.shape[1], 2)

        for vert in VERTICES.T:
            self.assertIn(vert, actual_vertices.T)

if __name__ == '__main__':
    unittest.main()
