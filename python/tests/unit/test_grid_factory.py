#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import the bempp module
import sys
sys.path.append("../..")
import bempp

import pytest

class TestGridFactory:
    def test_createStructuredGrid_throws_for_invalid_topology(self):
        pytest.raises((RuntimeError, ValueError),
            "bempp.GridFactory.createStructuredGrid('invalid', (0., 0.), (1., 2.), (4, 5))")

    def test_createStructuredGrid_throws_for_mismatched_dimensions_a(self):
        pytest.raises((RuntimeError, ValueError),
            "bempp.GridFactory.createStructuredGrid('triangular', (0., 0.), (1., 2., 3.), (4, 5))")

    def test_createStructuredGrid_throws_for_mismatched_dimensions_b(self):
        pytest.raises((RuntimeError, ValueError),
            "bempp.GridFactory.createStructuredGrid('triangular', (0., 0.), (1., 2.), (4, 5, 3))")

    def test_createStructuredGrid_throws_for_zero_elements_requested(self):
        pytest.raises((RuntimeError, ValueError),
            "bempp.GridFactory.createStructuredGrid('triangular', (0., 0.), (1., 2.), (0, 3))")

    def test_createGridFromConnectivityArrays_creates_correct_grid_without_domain_indices(self):
        import numpy as np
        vertices = np.array([[0, 0, 1, 1.2],
                             [0, 1, 0, 1.1],
                             [0, 0, 0, 0.5]])
        elementCorners = np.array([[0, 2],
                                   [1, 1],
                                   [2, 3],
                                   [-1, -1]])
        grid = bempp.GridFactory.createGridFromConnectivityArrays(
            'triangular', vertices, elementCorners)
        view = grid.leafView()
        v, ec, aux = view.getRawElementData()
        assert np.all(vertices == v)
        assert np.all(elementCorners == ec)
        assert aux.size == 0

    def test_createGridFromConnectivityArrays_creates_correct_grid_with_domain_indices(self):
        import numpy as np
        vertices = np.array([[0, 0, 1, 1.2],
                             [0, 1, 0, 1.1],
                             [0, 0, 0, 0.5]])
        elementCorners = np.array([[0, 2],
                                   [1, 1],
                                   [2, 3],
                                   [-1, -1]])
        domainIndices = [5, 7]
        grid = bempp.GridFactory.createGridFromConnectivityArrays(
            'triangular', vertices, elementCorners, domainIndices)
        view = grid.leafView()
        v, ec, aux, di = view.getRawElementData(returnDomainIndices=True)
        assert np.all(vertices == v)
        assert np.all(elementCorners == ec)
        assert aux.size == 0
        assert np.all(domainIndices == di)

    # TODO: test importGmshGrid -> create some small mesh for testing.
