#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import the bempp module
import sys
sys.path.append("..")
import bempp

import pytest

class TestEntity:
    def setup_method(self, method):
        self.grid = bempp.GridFactory.createStructuredGrid("triangular", (0., 0.), (1., 2.), (4, 5))
        self.view = self.grid.levelView(0)

    def get_entity(self, codim):
        it = self.view.entities(codim)
        it.next()
        return it.next()


    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_level_is_correct(self, codim):
        entity = self.get_entity(codim)
        assert entity.level() == 0


    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_type_returns_a_GeometryType_for_codim(self, codim):
        entity = self.get_entity(codim)
        assert isinstance(entity.type(), bempp.GeometryType)

    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_type_agrees_with_geometry_type_for_codim(self, codim):
        entity = self.get_entity(codim)
        assert entity.type() == entity.geometry().type()

    def test_type_is_correct_for_codim_0(self):
        entity = self.get_entity(0)
        assert entity.type() == bempp.GeometryType("simplex", 2)

    def test_type_is_correct_for_codim_1(self):
        entity = self.get_entity(1)
        assert entity.type() == bempp.GeometryType(1)

    def test_type_is_correct_for_codim_2(self):
        entity = self.get_entity(2)
        assert entity.type() == bempp.GeometryType(0)


    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_geometry_returns_a_Geometry_for_codim(self, codim):
        entity = self.get_entity(codim)
        assert isinstance(entity.geometry(), bempp.Geometry)


    @pytest.mark.parametrize("codim", (1, 2))
    def test_subEntityCount_does_not_exist_for_codim(self, codim):
        entity = self.get_entity(codim)
        assert not hasattr(entity, "subEntityCount")

    def test_subEntityCount_is_correct_for_codimSub_1(self):
        entity = self.get_entity(0)
        assert entity.subEntityCount(1) == 3   # three edges
    
    def test_subEntityCount_is_correct_for_codimSub_2(self):
        entity = self.get_entity(0)
        assert entity.subEntityCount(2) == 3   # three vertices

    def test_subEntityCount_is_0_for_codimSub_2(self):
        entity = self.get_entity(0)
        assert entity.subEntityCount(0) == 0

    def test_subEntityCount_is_0_for_codimSub_2(self):
        entity = self.get_entity(0)
        assert entity.subEntityCount(3) == 0

    @pytest.mark.parametrize("codimSub", (1, 2))
    def test_subEntityCount_agrees_with_number_of_iterations_for_codimSub(self, codimSub):
        entity = self.get_entity(0)
        allSubEntities = [e for e in entity.subEntities(codimSub)]
        assert entity.subEntityCount(codimSub) == len(allSubEntities)

    
    @pytest.mark.parametrize("codimSub", (0, 3))
    def test_subEntities_throws_for_codimSub(self, codimSub):
        entity = self.get_entity(0)
        pytest.raises(Exception, "entity.subEntities(codimSub)")

    def test_subEntities_returns_an_iterator_for_codimSub_1(self):
        entity = self.get_entity(0)
        assert isinstance(entity.subEntities(1), bempp.EntityIteratorCodim1)

    def test_subEntities_returns_an_iterator_for_codimSub_2(self):
        entity = self.get_entity(0)
        assert isinstance(entity.subEntities(2), bempp.EntityIteratorCodim2)


    def test_subEntities_next_returns_an_entity_for_codim_1(self):
        entity = self.get_entity(0)
        assert isinstance(entity.subEntities(1).next(), bempp.EntityCodim1)

    def test_subEntities_next_returns_an_entity_for_codim_2(self):
        entity = self.get_entity(0)
        assert isinstance(entity.subEntities(2).next(), bempp.EntityCodim2)

    @pytest.mark.parametrize("codimSub", (1, 2))
    def test_subEntities_next_parentGrid_is_correct(self, codimSub):
        entity = self.get_entity(0)
        subEntity = entity.subEntities(2).next()
        assert subEntity.parentGrid is self.grid

    
    def test_hasFather_is_correct(self):
        entity = self.get_entity(0)
        assert not entity.hasFather()

    def test_isLeaf_is_correct(self):
        entity = self.get_entity(0)
        assert entity.isLeaf()

    def test_isRegular_is_correct(self):
        entity = self.get_entity(0)
        assert entity.isRegular()


    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_parentGrid_is_correct(self, codim):
        entity = self.get_entity(codim)
        assert entity.parentGrid is self.grid

    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_entity_does_not_allow_to_delete_parent_grid(self, codim):
        entity = self.get_entity(codim)
        del self.grid
        del self.view
        entity.parentGrid.dim() # should not raise an exception

    # Functions related to refinement are not tested at present, since they are 
    # not yet implemented in Dune::FoamGrid.
