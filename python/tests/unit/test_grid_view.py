#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import the bempp module
import sys
sys.path.append("../..")
import bempp

import pytest

class BaseTestGridView:
    def test_indexSet_returns_an_IndexSet(self):
        assert isinstance(self.view.indexSet(), bempp.IndexSet)


    def test_entityCount_is_correct_for_codim_0(self):
        assert self.view.entityCount(0) == 4 * 5 * 2

    def test_entityCount_is_correct_for_codim_2(self):
        assert self.view.entityCount(2) == (4 + 1) * (5 + 1)

    def test_entityCount_is_zero_for_codim_3(self):
        assert self.view.entityCount(3) == 0


    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_entityCount_agrees_with_number_of_iterations_for_codim(self, codim):
        allEntities = [e for e in self.view.entities(codim)]
        assert self.view.entityCount(codim) == len(allEntities)

    def test_entityCount_is_correct_for_triangle(self):
        triangle = bempp.GeometryType("simplex", 2)
        assert self.view.entityCount(triangle) == 4 * 5 * 2


    def test_entities_returns_an_iterator_for_codim_0(self):
        assert isinstance(self.view.entities(0), bempp.EntityIteratorCodim0)

    def test_entities_returns_an_iterator_for_codim_1(self):
        assert isinstance(self.view.entities(1), bempp.EntityIteratorCodim1)

    def test_entities_returns_an_iterator_for_codim_2(self):
        assert isinstance(self.view.entities(2), bempp.EntityIteratorCodim2)

    def test_entities_throws_for_codim_3(self):
        pytest.raises(Exception, "self.view.entities(3)")

    def test_entities_next_returns_an_entity_for_codim_0(self):
        assert isinstance(self.view.entities(0).next(), bempp.EntityCodim0)

    def test_entities_next_returns_an_entity_for_codim_1(self):
        assert isinstance(self.view.entities(1).next(), bempp.EntityCodim1)

    def test_entities_next_returns_an_entity_for_codim_2(self):
        assert isinstance(self.view.entities(2).next(), bempp.EntityCodim2)


    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_containsEntity_returns_True_for_second_entity_of_codim(self, codim):
        it = self.view.entities(codim)
        e = it.next()
        e = it.next()
        assert self.view.containsEntity(e)


    def test_vtkWriter_with_no_arguments_returns_a_vtk_writer(self):
        assert isinstance(self.view.vtkWriter(), bempp.VtkWriter)
    
    def test_vtkWriter_returns_a_vtk_writer_for_dm_conforming(self):
        assert isinstance(self.view.vtkWriter("conforming"), bempp.VtkWriter)
    
    def test_vtkWriter_returns_a_vtk_writer_for_dm_nonconforming(self):
        assert isinstance(self.view.vtkWriter("nonconforming"), bempp.VtkWriter)

    def test_vtkWriter_throws_ValueError_for_invalid_dm(self):
        pytest.raises(ValueError, "self.view.vtkWriter('invalid')")


    def test_parentGrid_is_correct(self):
        assert self.view.parentGrid is self.grid

    def test_view_does_not_allow_to_delete_parent_grid(self):
        del self.grid
        self.view.parentGrid.dim() # should not raise an exception


class TestLevelGridView(BaseTestGridView):
    def setup_method(self, method):
        self.grid = bempp.GridFactory.createStructuredGrid("triangular", (0., 0.), (1., 2.), (4, 5))
        self.view = self.grid.levelView(0)

class TestLeafGridView(BaseTestGridView):
    def setup_method(self, method):
        self.grid = bempp.GridFactory.createStructuredGrid("triangular", (0., 0.), (1., 2.), (4, 5))
        self.view = self.grid.leafView()
