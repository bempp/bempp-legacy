import pytest
import numpy as np

class TestEigenConversion(object):

    m = 5
    n = 6

    def test_matrix_conversion_float32(self):

        from bempp.utils.eigen import _test_eigen_matrix_conversion_float32

        expected = np.random.rand(self.m,self.n).astype(dtype='float32')

        actual = _test_eigen_matrix_conversion_float32(expected)

        assert expected.dtype==actual.dtype
        assert expected.shape==actual.shape
        assert np.all(expected==actual)


    def test_matrix_conversion_float64(self):

        from bempp.utils.eigen import _test_eigen_matrix_conversion_float64

        expected = np.random.rand(self.m,self.n).astype(dtype='float64')

        actual = _test_eigen_matrix_conversion_float64(expected)

        assert expected.dtype==actual.dtype
        assert expected.shape==actual.shape
        assert np.all(expected==actual)


    def test_matrix_conversion_complex64(self):

        from bempp.utils.eigen import _test_eigen_matrix_conversion_complex64

        expected = (np.random.rand(self.m,self.n)+1j*np.random.rand(self.m,self.n)).astype(dtype='complex64')

        actual = _test_eigen_matrix_conversion_complex64(expected)

        assert expected.dtype==actual.dtype
        assert expected.shape==actual.shape
        assert np.all(expected==actual)

    def test_matrix_conversion_complex128(self):

        from bempp.utils.eigen import _test_eigen_matrix_conversion_complex128

        expected = (np.random.rand(self.m,self.n)+1j*np.random.rand(self.m,self.n)).astype(dtype='complex128')

        actual = _test_eigen_matrix_conversion_complex128(expected)

        assert expected.dtype==actual.dtype
        assert expected.shape==actual.shape
        assert np.all(expected==actual)

    def test_vector_conversion_float32(self):

        from bempp.utils.eigen import _test_eigen_vector_conversion_float32

        expected = np.random.rand(self.m).astype(dtype='float32')

        actual = _test_eigen_vector_conversion_float32(expected)

        assert expected.dtype==actual.dtype
        assert expected.shape==actual.shape
        assert np.all(expected==actual)


    def test_vector_conversion_float64(self):

        from bempp.utils.eigen import _test_eigen_vector_conversion_float64

        expected = np.random.rand(self.m).astype(dtype='float64')

        actual = _test_eigen_vector_conversion_float64(expected)

        assert expected.dtype==actual.dtype
        assert expected.shape==actual.shape
        assert np.all(expected==actual)


    def test_vector_conversion_complex64(self):

        from bempp.utils.eigen import _test_eigen_vector_conversion_complex64

        expected = (np.random.rand(self.m)+1j*np.random.rand(self.m)).astype(dtype='complex64')

        actual = _test_eigen_vector_conversion_complex64(expected)

        assert expected.dtype==actual.dtype
        assert expected.shape==actual.shape
        assert np.all(expected==actual)

    def test_vector_conversion_complex128(self):

        from bempp.utils.eigen import _test_eigen_vector_conversion_complex128

        expected = (np.random.rand(self.m)+1j*np.random.rand(self.m)).astype(dtype='complex128')

        actual = _test_eigen_vector_conversion_complex128(expected)

        assert expected.dtype==actual.dtype
        assert expected.shape==actual.shape
        assert np.all(expected==actual)
