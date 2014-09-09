from py.test import mark
from numpy import dtype


class TestBoundaryOperator(object):

    def test_instantiate_empty_boundary_operator(self):
        from bempp.assembly import BoundaryOperator
        BoundaryOperator("float32", "float32")

    @mark.parametrize('basis_type, result_type', [
        ("float32", "float32"),
        ("float32", "complex64"),
        ("float64", "complex128"),
        (dtype('float64'), dtype('complex128'))
    ])
    def test_dtypes(self, basis_type, result_type):
        from bempp.assembly import BoundaryOperator
        bop = BoundaryOperator(basis_type, result_type)
        assert bop.basis_type == basis_type
        assert bop.result_type == result_type

    @mark.parametrize('basis_type, result_type', [
        ("float64", "float32"),
        ("complex64", "float64"),
        ("float32", "complex128")
    ])
    def test_incompatible_dtypes(self, basis_type, result_type):
        from py.test import raises
        from bempp.assembly import BoundaryOperator
        with raises(TypeError):
            BoundaryOperator(basis_type, result_type)
