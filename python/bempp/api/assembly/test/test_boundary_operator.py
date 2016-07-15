from unittest import TestCase


class TestBoundaryOperator(TestCase):

    def setUp(self):
        import bempp.api

        grid = bempp.api.shapes.regular_sphere(2)
        self.domain = bempp.api.function_space(grid, "DP", 0)
        self.range_ = bempp.api.function_space(grid, "DP", 1)
        self.dual_to_range = bempp.api.function_space(grid, "DP", 2)
        parameters = bempp.api.common.global_parameters()
        parameters.assembly.use_super_spaces = False
        self._local_operator = bempp.api.operators.boundary.sparse.identity(
            self.domain, self.range_, self.dual_to_range)
        self._elementary_operator = bempp.api.operators.boundary.laplace.single_layer(
            self.domain, self.range_, self.dual_to_range, parameters=parameters)

    def test_elementary_boundary_operator_domain(self):
        self.assertTrue(self.domain.is_identical(
            self._elementary_operator.domain))

    def test_elementary_boundary_operator_range(self):
        self.assertTrue(self.range_.is_identical(
            self._elementary_operator.range))

    def test_elementary_boundary_operator_dual_to_range(self):
        self.assertTrue(self.dual_to_range.is_identical(
            self._elementary_operator.dual_to_range))

    def test_local_boundary_operator_domain(self):
        self.assertTrue(self.domain.is_identical(self._local_operator.domain))

    def test_local_boundary_operator_range(self):
        self.assertTrue(self.range_.is_identical(self._local_operator.range))

    def test_local_boundary_operator_dual_to_range(self):
        self.assertTrue(self.dual_to_range.is_identical(
            self._local_operator.dual_to_range))

    def test_sum_of_operators_is_sum_object(self):
        from bempp.api.assembly.boundary_operator import _SumBoundaryOperator

        self.assertIsInstance(self._local_operator + self._elementary_operator, _SumBoundaryOperator,
                              "Sum of operators should be of type _SumBoundaryOperator.")

    def test_product_of_operator_with_scalar_is_scaled_boundary_operator(self):
        from bempp.api.assembly.boundary_operator import _ScaledBoundaryOperator

        self.assertIsInstance(2.0 * self._elementary_operator, _ScaledBoundaryOperator,
                              "Sum of operators should be of type _ScaledBoundaryOperator.")

    def test_product_of_two_operators_is_product_operator(self):
        import bempp
        from bempp.api.assembly.boundary_operator import _ProductBoundaryOperator

        op = bempp.api.operators.boundary.laplace.single_layer(
            self.domain, self.domain, self.domain)

        self.assertIsInstance(op * op, _ProductBoundaryOperator,
                              "Product of two boundary operators should be _ProductBoundaryOperator.")

    def test_weak_form_of_local_operator_is_sparse_discrete_operator(self):
        from bempp.api.assembly.discrete_boundary_operator import SparseDiscreteBoundaryOperator

        self.assertIsInstance(self._local_operator.weak_form(),
                              SparseDiscreteBoundaryOperator)

    def test_weak_form_of_operator_sum_is_discrete_operator_sum(self):
        from bempp.api.assembly.discrete_boundary_operator import DiscreteBoundaryOperatorSum

        operator_sum = self._local_operator + self._elementary_operator

        self.assertIsInstance(operator_sum.weak_form(), DiscreteBoundaryOperatorSum,
                              "A DiscreteBoundaryOperatorSum instance is expected here.")

    def test_weak_form_of_scaled_operator_is_scaled_discrete_operator(self):
        from bempp.api.assembly.discrete_boundary_operator import ScaledDiscreteBoundaryOperator

        scaled_operator = 2.0 * self._elementary_operator
        weak_form = scaled_operator.weak_form()

        self.assertIsInstance(weak_form, ScaledDiscreteBoundaryOperator,
                              "A ScaledDiscreteBoundaryOperator instance is expected here. Actual type: " +
                              str(type(weak_form)))

    def test_weak_form_of_product_operator_is_product_discrete_operator(self):
        import bempp.api
        from bempp.api.assembly.discrete_boundary_operator import DiscreteBoundaryOperatorProduct

        op = bempp.api.operators.boundary.laplace.single_layer(
            self.domain, self.domain, self.domain)

        self.assertIsInstance((op * op).weak_form(), DiscreteBoundaryOperatorProduct,
                              "A DiscreteBoundaryOperatorProduct instance is expected.")


if __name__ == "__main__":
    from unittest import main

    main()
