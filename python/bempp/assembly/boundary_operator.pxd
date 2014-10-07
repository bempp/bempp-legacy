from libcpp cimport bool as cbool
from libcpp cimport complex
from libcpp.string cimport string
from bempp.space.space cimport SpaceVariants


cdef extern from "bempp/assembly/boundary_operator.hpp":
    cdef cppclass c_BoundaryOperator "Bempp::BoundaryOperator" [BASIS, RESULT]:
        c_BoundaryOperator()

cdef extern from "bempp/assembly/variants.hpp" namespace "Bempp":
    cdef cppclass BoundaryOpVariants:
        BoundaryOpVariants()
        void set[BASIS, RESULT](const c_BoundaryOperator[BASIS, RESULT] &_in) except+
        void set[BASIS, RESULT]() except+
        # Version to catches exceptions on assignemnt. Operator= is not
        # explicitly defined by bempp. 
        void assign "operator="(const BoundaryOpVariants&) except+
        string basisType() const
        string resultType() const
        SpaceVariants range() except +
        SpaceVariants dual_to_range "dualToRange"() except +
        SpaceVariants domain() except +
        BoundaryOpVariants operator+(const BoundaryOpVariants &_in) except+
        BoundaryOpVariants operator-(const BoundaryOpVariants &_in) except+
        BoundaryOpVariants operator*(const BoundaryOpVariants &_in) except+
        BoundaryOpVariants operator*(const float& _in) except+
        BoundaryOpVariants operator*(const double& _in) except+
        BoundaryOpVariants operator*(complex& _in) except+
        BoundaryOpVariants operator/(const float& _in) except+
        BoundaryOpVariants operator/(const double& _in) except+
        BoundaryOpVariants operator/(complex& _in) except+
        string label() const


cdef class BoundaryOperator:
    cdef BoundaryOpVariants impl_
