#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import the bempp module
import sys
sys.path.append("..")
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

    # TODO: test importGmshGrid -> create some small mesh for testing.