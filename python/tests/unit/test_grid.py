#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import the bempp module
import sys
sys.path.append("../..")
import bempp

import pytest

class TestGrid:
    def setup_method(self, method):
        self.grid = bempp.GridFactory.createStructuredGrid("triangular", (0., 0.), (1., 2.), (4, 5))

    def test_dim_of_triangular_grid_is_2(self):
        assert self.grid.dim() == 2

    def test_dimWorld_of_triangular_grid_is_3(self):
        assert self.grid.dimWorld() == 3

    def test_maxLevel_of_unrefined_grid_is_0(self):
        assert self.grid.maxLevel() == 0

    def test_levelView_returns_a_GridView_for_level_0(self):
        assert isinstance(self.grid.levelView(0), bempp.GridView)

    def test_levelView_of_unrefined_grid_throws_for_level_1(self):
        pytest.raises(Exception, "self.grid.levelView(1)")

    def test_leafView_returns_a_GridView(self):
        assert isinstance(self.grid.leafView(), bempp.GridView)
 
    def test_globalIdSet_returns_an_IdSet(self):
        assert isinstance(self.grid.globalIdSet(), bempp.IdSet)

