cdef extern from "bempp/grid/bempp_grid/test_grid.h" namespace "Bempp":
    cdef void c_test_grid "Bempp::test_grid"()

def test_grid():

    c_test_grid()



