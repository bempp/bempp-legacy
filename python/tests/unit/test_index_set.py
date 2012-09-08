#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import the bempp module
import sys
sys.path.append("../..")
import bempp

import pytest

class TestIndexSet:
    def setup_method(self, method):
        self.grid = bempp.GridFactory.createStructuredGrid("triangular", (0., 0.), (1., 2.), (4, 5))
        self.view = self.grid.leafView()
        self.indexSet = self.view.indexSet()

    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_entityIndex_returns_unique_values_for_codim(self, codim):
        allEntities = [e for e in self.view.entities(codim)]
        allIndices = [self.indexSet.entityIndex(e) for e in allEntities]
        assert len(allIndices) == len(set(allIndices))

    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_entityIndex_returns_correct_range_of_values_for_codim(self, codim):
        allEntities = [e for e in self.view.entities(codim)]
        allIndices = [self.indexSet.entityIndex(e) for e in allEntities]
        allIndices = set(allIndices)
        assert allIndices == set(range(len(allIndices)))

    @pytest.mark.parametrize("codimSub", (1, 2))
    def test_subEntityIndex_agrees_with_entityIndex_for_codimSub(self, codimSub):
	it = self.view.entities(0)
	it.next()
	entity = it.next()
	subIt = entity.subEntities(codimSub)
	subIt.next()
	subEntity = subIt.next()
	directIndex = self.indexSet.entityIndex(subEntity)
	indirectIndex = self.indexSet.subEntityIndex(entity, 1, codimSub)
        assert directIndex == indirectIndex

    def test_subEntityindex_returns_the_same_as_entityindex_for_codimSub_0(self):
	it = self.view.entities(0)
	it.next()
	entity = it.next()
	directIndex = self.indexSet.entityIndex(entity)
	indirectIndex = self.indexSet.subEntityIndex(entity, 0, 0)
        assert directIndex == indirectIndex

    def test_parentGridView_is_correct(self):
        assert self.indexSet.parentGridView is self.view