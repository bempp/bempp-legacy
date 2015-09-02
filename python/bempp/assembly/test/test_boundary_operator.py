from unittest import TestCase


class TestBoundaryOperator(TestCase):
    def setUp(self):
        import bempp

        grid = bempp.grid_from_sphere(2)
        self.domain = bempp.function_space(grid, "DP", 0)
        self.range_ = bempp.function_space(grid, "DP", 1)
        self.dual_to_range = bempp.function_space(grid, "DP", 2)
        self._local_operator = bempp.operators.boundary.sparse.identity(
            self.domain, self.range_, self.dual_to_range)
        self._elementary_operator = bempp.operators.boundary.laplace.single_layer(
            self.domain, self.range_, self.dual_to_range)

    def test_elementary_boundary_operator_domain(self):
        self.assertTrue(self.domain.is_identical(self._elementary_operator.domain))

    def test_elementary_boundary_operator_range(self):
        self.assertTrue(self.range_.is_identical(self._elementary_operator.range))

    def test_elementary_boundary_operator_dual_to_range(self):
        self.assertTrue(self.dual_to_range.is_identical(self._elementary_operator.dual_to_range))

    def test_local_boundary_operator_domain(self):
        self.assertTrue(self.domain.is_identical(self._local_operator.domain))

    def test_local_boundary_operator_range(self):
        self.assertTrue(self.range_.is_identical(self._local_operator.range))

    def test_local_boundary_operator_dual_to_range(self):
        self.assertTrue(self.dual_to_range.is_identical(self._local_operator.dual_to_range))

    def test_sum_of_operators_is_sum_object(self):
        from bempp.assembly.boundary_operator import _SumBoundaryOperator

        self.assertIsInstance(self._local_operator + self._elementary_operator, _SumBoundaryOperator,
                              "Sum of operators should be of type _SumBoundaryOperator.")

    def test_product_of_operator_with_scalar_is_scaled_boundary_operator(self):
        from bempp.assembly.boundary_operator import _ScaledBoundaryOperator

        self.assertIsInstance(2.0 * self._elementary_operator, _ScaledBoundaryOperator,
                              "Sum of operators should be of type _ScaledBoundaryOperator.")

    def test_product_of_two_operators_is_product_operator(self):
        import bempp
        from bempp.assembly.boundary_operator import _ProductBoundaryOperator

        op = bempp.operators.boundary.laplace.single_layer(self.domain, self.domain, self.domain)

        self.assertIsInstance(op * op, _ProductBoundaryOperator,
                              "Product of two boundary operators should be _ProductBoundaryOperator.")

    def test_weak_form_of_local_operator_is_sparse_discrete_operator(self):
        from bempp.assembly.discrete_boundary_operator import SparseDiscreteBoundaryOperator

        self.assertIsInstance(self._local_operator.weak_form(), SparseDiscreteBoundaryOperator)

    def test_weak_form_of_dense_operator_is_dense_discrete_operator(self):
        import bempp
        from bempp.assembly.discrete_boundary_operator import DenseDiscreteBoundaryOperator

        assembly_mode = bempp.global_parameters.assembly.boundary_operator_assembly_type
        bempp.global_parameters.assembly.boundary_operator_assembly_type = 'dense'
        self.assertIsInstance(self._elementary_operator.weak_form(), DenseDiscreteBoundaryOperator)
        bempp.global_parameters.assembly.boundary_operator_assembly_type = assembly_mode

    def test_weak_form_of_operator_sum_is_discrete_operator_sum(self):
        from bempp.utils.linear_operator import _SumLinearOperator

        operator_sum = self._local_operator + self._elementary_operator

        self.assertIsInstance(operator_sum.weak_form(), _SumLinearOperator,
                              "A _SumLinearOperator instance is expected here.")

    def test_weak_form_of_scaled_operator_is_scaled_discrete_operator(self):
        from bempp.utils.linear_operator import _ScaledLinearOperator

        scaled_operator = 2.0 * self._elementary_operator
        weak_form = scaled_operator.weak_form()

        self.assertIsInstance(weak_form, _ScaledLinearOperator,
                              "A _ScaledLinearOperator instance is expected here. Actual type: " +
                              str(type(weak_form)))

    def test_weak_form_of_product_operator_is_product_discrete_operator(self):
        import bempp
        from bempp.utils.linear_operator import _ProductLinearOperator

        op = bempp.operators.boundary.laplace.single_layer(self.domain, self.domain, self.domain)

        self.assertIsInstance((op * op).weak_form(), _ProductLinearOperator,
                              "A _ProductLinearOperator is expected.")


if __name__ == "__main__":
    from unittest import main

    main()
