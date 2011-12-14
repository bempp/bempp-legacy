#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import the bempp module
import sys
sys.path.append("..")
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

    def test_parentGrid_is_correct(self):
        assert self.idSet.parentGrid is self.grid    
