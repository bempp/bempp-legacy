<%
from bempp_operators import dtypes, compatible_dtypes, bops, ctypes
%>
from libcpp cimport bool as cbool
from libcpp.string cimport string
from bempp.utils cimport catch_exception
from bempp.space.space cimport SpaceVariants

# Declares complex type explicitly.
# Cython 0.20 will fail if templates are nested more than three-deep,
# as in shared_ptr[ c_Space[ complex[float] ] ]
cdef extern from "bempp/space/types.h":
% for ctype in dtypes.values():
%     if 'complex'  in ctype:
    ctypedef struct ${ctype}
%     endif
% endfor


cdef extern from "bempp/assembly/boundary_operator.hpp":
    cdef cppclass c_BoundaryOperator "Bempp::BoundaryOperator" [BASIS, RESULT]:
        c_BoundaryOperator()

cdef extern from "bempp/assembly/variants.hpp" namespace "Bempp":
    cdef cppclass BoundaryOpVariants:
        BoundaryOpVariants() except+
# Older cythons or cython with --gdb does not work well with templates
# So use the workaround below to explicitly define overloaded function with
# aliases to actual templates.
% for pybasis, cybasis in dtypes.items():
%     for pyresult, cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
        void set${cybasis}${cyresult} "set<${cybasis | ctypes}, ${cyresult | ctypes}>"(
                const c_BoundaryOperator[${cybasis}, ${cyresult}] &_in) except+
        void set${cybasis}${cyresult} "set<${cybasis | ctypes}, ${cyresult | ctypes}>"() except+
%         endif
%     endfor
% endfor
        # Version to catches exceptions on assignemnt. Operator= is not
        # explicitly defined by bempp.
        void assign "operator="(const BoundaryOpVariants&) except+
        string basisType() const
        string resultType() const
        cbool valid_range "validRange"()
        cbool valid_dual_to_range "validDualToRange"()
        cbool valid_domain "validDomain"()
        SpaceVariants range() except+catch_exception
        SpaceVariants dual_to_range "dualToRange"() except+catch_exception
        SpaceVariants domain() except+catch_exception
        BoundaryOpVariants operator+(const BoundaryOpVariants &_in) except+catch_exception
        BoundaryOpVariants operator-(const BoundaryOpVariants &_in) except+catch_exception
        BoundaryOpVariants operator*(const BoundaryOpVariants &_in) except+catch_exception
        BoundaryOpVariants operator*(const float& _in) except +catch_exception
        BoundaryOpVariants operator*(const double& _in) except +catch_exception
        BoundaryOpVariants operator*(const float complex& _in) except +catch_exception
        BoundaryOpVariants operator*(const double complex& _in) except +catch_exception
        BoundaryOpVariants operator/(const float& _in) except +catch_exception
        BoundaryOpVariants operator/(const double& _in) except +catch_exception
        BoundaryOpVariants operator/(const float complex& _in) except +catch_exception
        BoundaryOpVariants operator/(const double complex& _in) except +catch_exception
        string label() const


cdef class BoundaryOperator:
    cdef BoundaryOpVariants impl_
