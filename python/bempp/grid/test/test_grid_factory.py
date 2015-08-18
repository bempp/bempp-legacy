"""Unit tests for the GridFactory."""

import unittest
import numpy as np


class TestGridFactory(unittest.TestCase):
    """Unit Tests for the GridFactory."""

    def setUp(self):
        from bempp import grid_from_sphere
        from bempp import GridFactory

        sphere = grid_from_sphere(3)

        self.vertices = sphere.leaf_view.vertices
        self.elements = sphere.leaf_view.elements

        self.n_elements = sphere.leaf_view.entity_count(0)
        self.n_vertices = sphere.leaf_view.entity_count(2)

        self.elem_permutation = range(self.n_elements)[::-1]
        self.vertex_permutation = range(self.n_vertices)[::-1]

        self.domain_indices = np.random.randint(10, size=self.n_elements)

        factory = GridFactory()
        for index in self.vertex_permutation:
            factory.insert_vertex(self.vertices[:, index])
        for i, index in enumerate(self.elem_permutation):
            elem = [self.vertex_permutation[self.elements[0, index]],
                    self.vertex_permutation[self.elements[1, index]],
                    self.vertex_permutation[self.elements[2, index]]]
            factory.insert_element(elem, self.domain_indices[i])

        self.grid = factory.finalize()

    def test_grid_size(self):

        self.assertEqual(self.grid.leaf_view.entity_count(0), self.n_elements)
        self.assertEqual(self.grid.leaf_view.entity_count(2), self.n_vertices)

    def test_element_insertion_indices(self):

        for elem in self.grid.leaf_view.entity_iterator(0):
            insertion_index = self.grid.element_insertion_index(elem)
            original_element_index = self.n_elements - insertion_index - 1
            corners = elem.geometry.corners
            for i in range(3):
                original_vertex_index = self.elements[i, original_element_index]
                self.assertAlmostEqual(np.max(np.abs(corners[:, i] - self.vertices[:, original_vertex_index])), 0, 13)

    def test_vertex_insertion_indices(self):

        for vertex in self.grid.leaf_view.entity_iterator(2):
            insertion_index = self.grid.vertex_insertion_index(vertex)
            original_index = self.n_vertices - insertion_index - 1
            self.assertAlmostEqual(np.max(np.abs(vertex.geometry.corners[:, 0] - self.vertices[:, original_index])), 0,
                                   13)

    def test_domain_indices(self):

        for elem in self.grid.leaf_view.entity_iterator(0):
            insertion_index = self.grid.element_insertion_index(elem)
            self.assertEqual(elem.domain, self.domain_indices[insertion_index])


if __name__ == '__main__':
    unittest.main()
