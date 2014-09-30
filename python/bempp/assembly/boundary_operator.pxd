from libcpp cimport bool as cbool
from libcpp.string cimport string
from bempp.space.space cimport SpaceVariants


cdef extern from "bempp/assembly/boundary_operator.hpp":
    cdef cppclass c_BoundaryOperator "Bempp::BoundaryOperator" [BASIS, RESULT]:
        c_BoundaryOperator()

cdef extern from "bempp/assembly/variants.hpp" namespace "Bempp":
    cdef cppclass BoundaryOpVariants:
        BoundaryOpVariants()
        void set[BASIS, RESULT](const c_BoundaryOperator[BASIS, RESULT] &_in)
        void set[BASIS, RESULT]()
        string basisType() const
        string resultType() const
        SpaceVariants range() except +
        SpaceVariants dual_to_range "dualToRange"() except +
        SpaceVariants domain() except +

cdef class BoundaryOperator:
    cdef BoundaryOpVariants impl_
