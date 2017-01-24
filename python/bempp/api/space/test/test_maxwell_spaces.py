"""Test Maxwell spaces."""

from unittest import TestCase
import unittest
import bempp.api
import numpy as np

#pylint: disable=invalid-name
#pylint: disable=protected-access
#pylint: disable=too-many-locals

class TestMaxwellSpaces(TestCase):
    """Test Maxwell spaces."""

    requiresgmsh = unittest.skipIf(
        bempp.api.GMSH_PATH is None, reason="Needs GMSH")

    def setUp(self):

        self._grid = bempp.api.shapes.sphere(h=0.2)

    @requiresgmsh
    def test_rwg_space_is_correctly_scaled(self):
        """RWG space is correctly scaled."""

        space_rt = bempp.api.function_space(self._grid, "RT", 0)
        space_rwg = bempp.api.function_space(self._grid, "RWG", 0)

        for element in self._grid.leaf_view.entity_iterator(0):
            dofs_rt, weights_rt = space_rt.get_global_dofs(
                element, dof_weights=True)
            dofs_rwg, weights_rwg = space_rwg.get_global_dofs(
                element, dof_weights=True)
            self.assertEqual(len(weights_rt), len(weights_rwg))
            edge_volumes = np.array(
                [edge.geometry.volume
                 for edge in element.sub_entity_iterator(1)])
            self.assertEqual(dofs_rt, dofs_rwg)
            for ind in range(3):
                self.assertAlmostEqual(
                    weights_rt[ind] * edge_volumes[ind], weights_rwg[ind])

    @requiresgmsh
    def test_maxwell_barycentric_rwg_mass_agrees_with_rwg_mass(self):
        """Barycentric RWG space is correct."""
        space_rwg = bempp.api.function_space(self._grid, "RWG", 0)
        space_rwg_b = bempp.api.function_space(self._grid, "B-RWG", 0)

        ident_rwg = bempp.api.operators.boundary.sparse.identity(
            space_rwg, space_rwg, space_rwg)
        ident_rwg_b = bempp.api.operators.boundary.sparse.identity(
            space_rwg_b, space_rwg_b, space_rwg_b)

        mat_ident_rwg = ident_rwg.weak_form().sparse_operator
        mat_ident_rwg_b = ident_rwg_b.weak_form().sparse_operator
        max_diff = np.max(np.abs((mat_ident_rwg - mat_ident_rwg_b).data))
        self.assertAlmostEqual(max_diff, 0)

    @requiresgmsh
    def test_nedelec_barycentric_space_(self):
        """Test Nedelec barycentric space."""
        nc_space = bempp.api.function_space(self._grid, "B-NC", 0)
        rt_space = bempp.api.function_space(self._grid, "B-RT", 0)

        op1 = bempp.api.operators.boundary.sparse.identity(
            nc_space, nc_space, nc_space)
        op2 = bempp.api.operators.boundary.sparse._maxwell_identity(
            rt_space, nc_space, nc_space)

        m1 = bempp.api.as_matrix(op1.weak_form())
        m2 = bempp.api.as_matrix(op2.weak_form())
        max_diff = np.max(np.abs((m1 + m2).data))
        self.assertAlmostEqual(max_diff, 0)

    @requiresgmsh
    def test_nedelec_barycentric_space_agrees_with_non_barycentric_space(self):
        """Nedelec barycentric space agrees with non-barycentric space."""

        nc_space = bempp.api.function_space(self._grid, "NC", 0)
        nc_space_bary = bempp.api.function_space(self._grid, "B-NC", 0)

        op1 = bempp.api.operators.boundary.sparse.identity(
            nc_space, nc_space, nc_space)
        op2 = bempp.api.operators.boundary.sparse.identity(
            nc_space_bary, nc_space_bary, nc_space_bary)

        m1 = bempp.api.as_matrix(op1.weak_form())
        m2 = bempp.api.as_matrix(op2.weak_form())
        max_diff = np.max(np.abs((m1 - m2).data))
        self.assertAlmostEqual(max_diff, 0)

    @requiresgmsh
    def test_scaled_nedelec_barycentric_spaceself(self):
        """Test scaled Nedelec barycentric space."""
        nc_space = bempp.api.function_space(self._grid, "SNC", 0)
        nc_space_bary = bempp.api.function_space(self._grid, "B-SNC", 0)

        op1 = bempp.api.operators.boundary.sparse.identity(
            nc_space, nc_space, nc_space)
        op2 = bempp.api.operators.boundary.sparse.identity(
            nc_space_bary, nc_space_bary, nc_space_bary)

        m1 = bempp.api.as_matrix(op1.weak_form())
        m2 = bempp.api.as_matrix(op2.weak_form())
        max_diff = np.max(np.abs((m1 - m2).data))
        self.assertAlmostEqual(max_diff, 0)

    @requiresgmsh
    def test_nedelec_space_is_tangential_conforming(self):
        """Nedelec space is tangential conforming."""
        grid = self._grid
        space = bempp.api.function_space(grid, "NC", 0)

        index_set = grid.leaf_view.index_set()
        #pylint: disable=unnecessary-lambda
        elements = sorted(list(grid.leaf_view.entity_iterator(0)),
                          key=lambda e: index_set.entity_index(e))

        for dof_index in range(space.global_dof_count):

            local_dofs = space.global_to_local_dofs([dof_index])
            local_dof_0 = local_dofs[0][0][0]
            local_dof_1 = local_dofs[0][0][1]
            mult_0 = local_dofs[1][0][0]
            mult_1 = local_dofs[1][0][1]

            elem_0 = elements[local_dof_0.entity_index]
            elem_1 = elements[local_dof_1.entity_index]

            # Compute tangents and normals

            local_dof_index_0 = local_dof_0.dof_index
            local_dof_index_1 = local_dof_1.dof_index

            if local_dof_index_0 == 0:
                p_0 = np.array([[0.5], [0]])
            if local_dof_index_0 == 1:
                p_0 = np.array([[0], [0.5]])
            if local_dof_index_0 == 2:
                p_0 = np.array([[0.5], [0.5]])

            if local_dof_index_1 == 0:
                p_1 = np.array([[0.5], [0]])
            if local_dof_index_1 == 1:
                p_1 = np.array([[0], [0.5]])
            if local_dof_index_1 == 2:
                p_1 = np.array([[0.5], [0.5]])

            edge_0 = list(elem_0.sub_entity_iterator(1))[local_dof_0.dof_index]
            tangent_0 = (edge_0.geometry.corners[:, 1] -
                         edge_0.geometry.corners[:, 0])

            normal_0 = np.cross(tangent_0, elem_0.geometry.normals(p_0)[:, 0])
            normal_1 = np.cross(tangent_0, elem_1.geometry.normals(p_1)[:, 0])
            normal_0 = normal_0 / np.linalg.norm(normal_0)
            normal_1 = normal_1 / np.linalg.norm(normal_1)

            shapeset_0 = space.shapeset(elem_0)
            shapeset_1 = space.shapeset(elem_1)

            val_0 = mult_0 * np.squeeze(
                shapeset_0.evaluate(p_0, local_dof_index_0, values=True))
            val_1 = mult_1 * np.squeeze(
                shapeset_1.evaluate(p_1, local_dof_index_1, values=True))

            inverse_jacob_transpose_0 = \
                elem_0.geometry.jacobian_inverses_transposed(p_0)[0]
            inverse_jacob_transpose_1 = \
                elem_1.geometry.jacobian_inverses_transposed(p_1)[0]

            elem_val_0 = inverse_jacob_transpose_0.dot(val_0)
            elem_val_1 = inverse_jacob_transpose_1.dot(val_1)

            elem_tangent_0 = elem_val_0.dot(tangent_0)
            elem_tangent_1 = elem_val_1.dot(tangent_0)

            self.assertAlmostEqual(elem_tangent_0, elem_tangent_1)

    @requiresgmsh
    def test_rt_space_is_normal_conforming(self):
        """RT Space is normal conforming."""
        grid = self._grid
        space = bempp.api.function_space(grid, "RT", 0)

        index_set = grid.leaf_view.index_set()
        #pylint: disable=unnecessary-lambda
        elements = sorted(list(grid.leaf_view.entity_iterator(0)),
                          key=lambda e: index_set.entity_index(e))

        for dof_index in range(space.global_dof_count):

            local_dofs = space.global_to_local_dofs([dof_index])
            local_dof_0 = local_dofs[0][0][0]
            local_dof_1 = local_dofs[0][0][1]
            mult_0 = local_dofs[1][0][0]
            mult_1 = local_dofs[1][0][1]

            elem_0 = elements[local_dof_0.entity_index]
            elem_1 = elements[local_dof_1.entity_index]

            # Compute tangents and normals

            local_dof_index_0 = local_dof_0.dof_index
            local_dof_index_1 = local_dof_1.dof_index

            if local_dof_index_0 == 0:
                p_0 = np.array([[0.5], [0]])
            if local_dof_index_0 == 1:
                p_0 = np.array([[0], [0.5]])
            if local_dof_index_0 == 2:
                p_0 = np.array([[0.5], [0.5]])

            if local_dof_index_1 == 0:
                p_1 = np.array([[0.5], [0]])
            if local_dof_index_1 == 1:
                p_1 = np.array([[0], [0.5]])
            if local_dof_index_1 == 2:
                p_1 = np.array([[0.5], [0.5]])

            edge_0 = list(elem_0.sub_entity_iterator(1))[local_dof_0.dof_index]
            edge_1 = list(elem_1.sub_entity_iterator(1))[local_dof_1.dof_index]
            tangent_0 = (edge_0.geometry.corners[:, 1] -
                         edge_0.geometry.corners[:, 0])
            tangent_1 = (edge_1.geometry.corners[:, 1] -
                         edge_1.geometry.corners[:, 0])

            normal_0 = np.cross(tangent_0, elem_0.geometry.normals(p_0)[:, 0])
            normal_1 = np.cross(tangent_1, elem_1.geometry.normals(p_1)[:, 0])
            normal_0 = normal_0 / np.linalg.norm(normal_0)
            normal_1 = normal_1 / np.linalg.norm(normal_1)

            shapeset_0 = space.shapeset(elem_0)
            shapeset_1 = space.shapeset(elem_1)

            val_0 = mult_0 * np.squeeze(
                shapeset_0.evaluate(p_0, local_dof_index_0, values=True))
            val_1 = mult_1 * np.squeeze(
                shapeset_1.evaluate(p_1, local_dof_index_1, values=True))

            jacob_transpose_0 = elem_0.geometry.jacobians_transposed(p_0)[0]
            jacob_transpose_1 = elem_1.geometry.jacobians_transposed(p_1)[0]

            elem_val_0 = (jacob_transpose_0.T.dot(val_0) /
                          elem_0.geometry.integration_elements(p_0)[0])
            elem_val_1 = (jacob_transpose_1.T.dot(val_1) /
                          elem_1.geometry.integration_elements(p_1)[0])

            elem_normal_0 = elem_val_0.dot(normal_0)
            elem_normal_1 = elem_val_1.dot(normal_1)

            self.assertAlmostEqual(elem_normal_0, elem_normal_1)


    #requiresgmsh
    def test_twisted_rt_space_agrees_with_nedelec_space(self):
        """Test that Nedelec spaces and twisted RT spaces agree."""

        grid = bempp.api.shapes.cube(h=.1)
        space_rt = bempp.api.function_space(grid, "RT", 0)
        space_nc = bempp.api.function_space(grid, "NC", 0)

        op1 = bempp.api.operators.boundary.sparse._maxwell_identity(
            space_rt, space_rt, space_nc)
        op2 = bempp.api.operators.boundary.sparse.identity(
            space_nc, space_rt, space_nc)

        m1 = bempp.api.as_matrix(op1.weak_form())
        m2 = bempp.api.as_matrix(op2.weak_form())
        diff = m1 + m2

        self.assertAlmostEqual(np.max(np.abs(diff.data)), 0)

    #requiresgmsh
    def test_twisted_rwg_space_agrees_with_scaled_nedelec_space(self):
        """Test that scaled Nedelec spaces and twisted RWG spaces agree."""

        grid = bempp.api.shapes.cube(h=.1)
        space_rwg = bempp.api.function_space(grid, "RWG", 0)
        space_nc = bempp.api.function_space(grid, "SNC", 0)

        op1 = bempp.api.operators.boundary.sparse._maxwell_identity(
            space_rwg, space_rwg, space_nc)
        op2 = bempp.api.operators.boundary.sparse.identity(
            space_nc, space_rwg, space_nc)

        m1 = bempp.api.as_matrix(op1.weak_form())
        m2 = bempp.api.as_matrix(op2.weak_form())

        # Need that they are negative to each other
        diff = m1 + m2

        self.assertAlmostEqual(np.max(np.abs(diff.data)), 0)



if __name__ == "__main__":
    from unittest import main

    main()
