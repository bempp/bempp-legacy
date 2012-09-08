#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import the bempp module
import sys
sys.path.append("../..")
import bempp

import pytest

class TestIdSet:
    def setup_method(self, method):
        self.grid = bempp.GridFactory.createStructuredGrid("triangular", (0., 0.), (1., 2.), (4, 5))
        self.view = self.grid.leafView()
        self.idSet = self.grid.globalIdSet()

    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_entityId_returns_unique_values_for_codim(self, codim):
        allEntities = [e for e in self.view.entities(codim)]
        allIds = [self.idSet.entityId(e) for e in allEntities]
        assert len(allIds) == len(set(allIds))

    @pytest.mark.parametrize("codimSub", (1, 2))
    def test_subEntityId_agrees_with_entityId_for_codimSub(self, codimSub):
	it = self.view.entities(0)
	it.next()
	entity = it.next()
	subIt = entity.subEntities(codimSub)
	subIt.next()
	subEntity = subIt.next()
	directId = self.idSet.entityId(subEntity)
	indirectId = self.idSet.subEntityId(entity, 1, codimSub)
        assert directId == indirectId

    def test_subEntityId_returns_the_same_as_entityId_for_codimSub_0(self):
	it = self.view.entities(0)
	it.next()
	entity = it.next()
	directId = self.idSet.entityId(entity)
	indirectId = self.idSet.subEntityId(entity, 0, 0)
        assert directId == indirectId

    def test_parentGrid_is_correct(self):
        assert self.idSet.parentGrid is self.grid    
