<%
from bempp_operators import dtypes, compatible_dtypes, bops, ctypes
%>
from libcpp cimport bool as cbool
from libcpp.string cimport string
from bempp.utils cimport catch_exception, complex_float, complex_double
from bempp.space.space cimport SpaceVariants
from bempp.assembly.discrete_boundary_operator cimport c_DiscreteBoundaryOperator
from bempp.assembly.discrete_boundary_operator cimport DiscreteBoundaryOperator
from bempp.assembly.discrete_boundary_operator cimport DiscreteBoundaryOperatorBase
from bempp.utils cimport shared_ptr
from bempp.space.space cimport Space
from bempp.utils cimport ParameterList
from libcpp cimport bool as cbool

cdef extern from "bempp/assembly/boundary_operator.hpp":
    cdef cppclass c_BoundaryOperator "Bempp::BoundaryOperator" [BASIS, RESULT]:
        c_BoundaryOperator()

cdef extern from "bempp/assembly/py_boundary_operator_variants.hpp" namespace "Bempp":
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
        string label() const

    cdef shared_ptr[c_DiscreteBoundaryOperator[ResultType]] _boundary_operator_variant_weak_form "Bempp::boundary_op_variant_weak_form" [BasisFunctionType,ResultType] (const BoundaryOpVariants& variant)

    cdef c_BoundaryOperator[BASIS, RESULT] _py_get_boundary_operator[BASIS,RESULT](const BoundaryOpVariants&);

cdef extern from "bempp/assembly/py_discrete_operator_support.hpp" namespace "Bempp":
    cdef object py_get_sparse_from_discrete_operator[VALUE](shared_ptr[c_DiscreteBoundaryOperator[VALUE]])

cdef class BoundaryOperatorBase:
    cdef object _basis_type
    cdef object _result_type

cdef class GeneralBoundaryOperator(BoundaryOperatorBase):
    cdef BoundaryOpVariants impl_
    cdef ParameterList _parameters
    cdef cbool _is_sparse

cdef class DenseBoundaryOperator(GeneralBoundaryOperator):
    pass

cdef class SparseBoundaryOperator(GeneralBoundaryOperator):
    cdef Space _domain
    cdef Space _range
    cdef Space _dual_to_range

