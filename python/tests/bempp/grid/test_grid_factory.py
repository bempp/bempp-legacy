import pytest
import numpy as np


class TestGridFactory:

    from bempp import grid_from_sphere

    grid = grid_from_sphere(3)

    vertices = grid.leaf_view.vertices
    elements = grid.leaf_view.elements

    n_elements = grid.leaf_view.entity_count(0)
    n_vertices = grid.leaf_view.entity_count(2)

    elem_permutation = range(n_elements)[::-1]
    vertex_permutation = range(n_vertices)[::-1]
    
    domain_indices = np.random.randint(10,size=n_elements)

    @pytest.fixture
    def grid(self):

        from bempp import GridFactory


        factory = GridFactory()
        for index in self.vertex_permutation:
            factory.insert_vertex(self.vertices[:,index])
        for i,index in enumerate(self.elem_permutation):
            elem = [self.vertex_permutation[self.elements[0,index]],
                    self.vertex_permutation[self.elements[1,index]],
                    self.vertex_permutation[self.elements[2,index]]]
            factory.insert_element(elem,self.domain_indices[i])

        return factory.finalize()

    def test_grid_size(self,grid):

        assert grid.leaf_view.entity_count(0) == self.n_elements
        assert grid.leaf_view.entity_count(2) == self.n_vertices

    def test_element_insertion_indices(self,grid):

        for elem in grid.leaf_view.entity_iterator(0):
            insertion_index = grid.element_insertion_index(elem)
            original_element_index = self.n_elements-insertion_index-1
            corners = elem.geometry.corners
            for i in range(3):
                original_vertex_index = self.elements[i,original_element_index]
                assert np.max(np.abs(corners[:,i]-self.vertices[:,original_vertex_index]))<1E-15

    def test_vertex_insertion_indices(self,grid):
        
        for vertex in grid.leaf_view.entity_iterator(2):
            insertion_index = grid.vertex_insertion_index(vertex)
            original_index = self.n_vertices-insertion_index-1
            assert np.max(np.abs(vertex.geometry.corners[:,0]-self.vertices[:,original_index]))<1E-15

    def test_domain_indices(self,grid):

        for elem in grid.leaf_view.entity_iterator(0):
            insertion_index = grid.element_insertion_index(elem)
            assert elem.domain == self.domain_indices[insertion_index]









