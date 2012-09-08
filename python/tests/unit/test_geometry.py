#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import the bempp module
import sys
sys.path.append("../..")
import bempp
import numpy as np

import pytest

class TestGeometry:
    def setup_method(self, method):
        self.grid = bempp.GridFactory.createStructuredGrid("triangular", (0., 0.), (1., 2.), (4, 5))
        self.view = self.grid.levelView(0)

    def get_entity(self, codim):
        it = self.view.entities(codim)
        it.next()
        return it.next()

    def get_local(self, codim, n_points):
        d = np.arange(self.grid.dim() - codim)[:,np.newaxis] # row index
        p = np.arange(n_points)[np.newaxis,:] # column index
        local = 0.1 * (d + 1) + 0.01 * (p + 1)
        return local

    def get_global(self, n_points):
        d = np.arange(self.grid.dimWorld())[:,np.newaxis] # row index
        p = np.arange(n_points)[np.newaxis,:] # column index
        global_ = 0.1 * (d + 1) + 0.01 * (p + 1)
        return global_


    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_affine_is_correct_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        assert geometry.affine()


    def test_cornerCount_is_correct_for_codim_0(self):
        geometry = self.get_entity(0).geometry()
        assert geometry.cornerCount() == 3

    def test_cornerCount_is_correct_for_codim_1(self):
        geometry = self.get_entity(1).geometry()
        assert geometry.cornerCount() == 2

    def test_cornerCount_is_correct_for_codim_2(self):
        geometry = self.get_entity(2).geometry()
        assert geometry.cornerCount() == 1


    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_corners_is_array_of_proper_shape_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        corners = geometry.corners()
        # We join these tests in a single function because the second one 
        # doesn't make sense anyway if the first one fails.
        assert isinstance(corners, np.ndarray)
        assert corners.shape == (self.grid.dimWorld(), geometry.cornerCount())

    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_corners_is_finite_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        corners = geometry.corners()
        assert np.all(np.isfinite(corners))


    @pytest.mark.parametrize("codim", (0, 1))
    def test_local2global_is_array_of_proper_shape_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        local = self.get_local(codim, n_points=5)
        global_ = geometry.local2global(local)
        assert isinstance(global_, np.ndarray)
        assert global_.shape == (self.grid.dimWorld(), 5)

    @pytest.mark.parametrize("codim", (0, 1))
    def test_local2global_is_finite_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        local = self.get_local(codim, n_points=5)
        global_ = geometry.local2global(local)
        assert np.all(np.isfinite(global_))

    @pytest.mark.parametrize("codim", (0, 1))
    def test_local2global_is_identical_for_inputs_being_a_2d_array_and_a_nested_list_and_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        local = self.get_local(codim, n_points=5)
        global_array = geometry.local2global(local)
        global_list = geometry.local2global(local.tolist())
        assert np.all(global_array == global_list)

    @pytest.mark.parametrize("codim", (0, 1))
    def test_local2global_throws_for_input_of_incorrect_dimensions_and_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        local = np.zeros((5, 4))
        pytest.raises((RuntimeError, ValueError), "geometry.local2global(local)")

    @pytest.mark.parametrize("codim", (0, 1))
    def test_local2global_throws_for_input_of_too_large_a_rank_and_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        local = np.zeros((5, 4, 1))
        pytest.raises((RuntimeError, ValueError, TypeError), "geometry.local2global(local)")

    @pytest.mark.parametrize("codim", (0, 1))
    def test_local2global_does_not_throw_for_zero_input_points_and_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        local = self.get_local(codim, n_points=0)
        global_ = geometry.local2global(local) # should not throw


    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_global2local_is_array_of_proper_shape_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        global_ = self.get_global(n_points=5)
        local = geometry.global2local(global_)
        assert isinstance(local, np.ndarray)
        assert local.shape == (self.grid.dim() - codim, 5)

    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_global2local_is_finite_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        global_ = self.get_global(n_points=5)
        local = geometry.global2local(global_)
        assert np.all(np.isfinite(local))


    @pytest.mark.parametrize("codim", (0, 1))
    def test_integrationElement_is_array_of_proper_shape_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        local = self.get_local(codim, n_points=5)
        integrationElements = geometry.integrationElements(local)
        assert isinstance(integrationElements, np.ndarray)
        assert integrationElements.shape == (5,)

    @pytest.mark.parametrize("codim", (0, 1))
    def test_integrationElement_is_finite_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        local = self.get_local(codim, n_points=5)
        integrationElements = geometry.integrationElements(local)
        assert np.all(np.isfinite(integrationElements))


    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_volume_is_nonnegative_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        volume = geometry.volume()
        assert volume >= 0


    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_center_is_array_of_proper_shape_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        center = geometry.center()
        # We join these tests in a single function because the second one 
        # doesn't make sense anyway if the first one fails.
        assert isinstance(center, np.ndarray)
        assert center.shape == (self.grid.dimWorld(),)

    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_center_is_finite_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        center = geometry.center()
        assert np.all(np.isfinite(center))


    @pytest.mark.parametrize("codim", (0, 1))
    def test_jacobianTransposed_is_array_of_proper_shape_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        local = self.get_local(codim, n_points=5)
        jacobiansTransposed = geometry.jacobiansTransposed(local)
        assert isinstance(jacobiansTransposed, np.ndarray)
        assert jacobiansTransposed.shape == \
            (self.grid.dim() - codim, self.grid.dimWorld(), 5)

    @pytest.mark.parametrize("codim", (0, 1))
    def test_jacobianTransposed_is_finite_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        local = self.get_local(codim, n_points=5)
        jacobiansTransposed = geometry.jacobiansTransposed(local)
        assert np.all(np.isfinite(jacobiansTransposed))


    @pytest.mark.parametrize("codim", (0, 1))
    def test_jacobianInverseTransposed_is_array_of_proper_shape_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        local = self.get_local(codim, n_points=5)
        jacobianInversesTransposed = geometry.jacobianInversesTransposed(local)
        assert isinstance(jacobianInversesTransposed, np.ndarray)
        assert jacobianInversesTransposed.shape == \
            (self.grid.dimWorld(), self.grid.dim() - codim, 5)

    @pytest.mark.parametrize("codim", (0, 1))
    def test_jacobianInverseTransposed_is_finite_for_codim(self, codim):
        geometry = self.get_entity(codim).geometry()
        local = self.get_local(codim, n_points=5)
        jacobianInversesTransposed = geometry.jacobianInversesTransposed(local)
        assert np.all(np.isfinite(jacobianInversesTransposed))


    @pytest.mark.parametrize("codim", (0, 1, 2))
    def test_parentEntity_is_correct(self, codim):
        entity = self.get_entity(codim)
        geometry = entity.geometry()
        assert geometry.parentEntity is entity