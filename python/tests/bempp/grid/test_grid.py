from py.test import fixture, mark
import numpy as np

class TestGrid(object):

    vertices = np.array([[0,1,1,0],
                         [0,0,1,1],
                         [0,0,0,0]])
    elements = np.array([[0,1],
                         [1,2],
                         [3,3]])

    def test_structured_grid(self):
        from bempp import structured_grid

        nx = 4
        ny = 5

        grid =  structured_grid(
            (0., 0.),
            ('1.', '2.'),
            (nx, ny)
        )

        box = grid.bounding_box

        assert np.all(box==np.array([[0,0,0],[1,2,0]]))
        assert grid.leaf_view.entity_count(0)==2*nx*ny

    def test_grid_from_element_data(self):

        from bempp import grid_from_element_data

        grid = grid_from_element_data(self.vertices,self.elements)
        
        box = grid.bounding_box
        assert np.all(box==np.array([[0,0,0],[1,1,0]]))
        assert grid.leaf_view.entity_count(0)==2
        assert grid.leaf_view.entity_count(2)==4
        
        actual_vertices = grid.leaf_view.vertices
        actual_elements = grid.leaf_view.elements

        assert actual_vertices.shape[1]==4
        assert actual_elements.shape[1]==2

        for vert in self.vertices.T:
            assert vert in actual_vertices.T



        

