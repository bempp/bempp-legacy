

cdef extern from "bempp/cuda/test_cuda.hpp" namespace "Bempp":
    cdef void c_test_cuda "Bempp::test_cuda"()

def test_cuda():

    c_test_cuda()

