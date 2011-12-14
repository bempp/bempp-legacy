#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import the bempp module
import sys
sys.path.append("..")
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

    def test_parentGridView_is_correct(self):
        assert self.indexSet.parentGridView is self.view