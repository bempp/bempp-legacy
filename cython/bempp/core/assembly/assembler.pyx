from bempp.core.utils cimport shared_ptr
from bempp.core.utils cimport complex_double
from bempp.core.utils cimport catch_exception
from bempp.core.utils cimport c_ParameterList, ParameterList
from bempp.core.space.space cimport c_Space, Space
from .discrete_boundary_operator cimport c_DiscreteBoundaryOperator
from .discrete_boundary_operator cimport RealDiscreteBoundaryOperator
from .discrete_boundary_operator cimport ComplexDiscreteBoundaryOperator
from cython.operator cimport dereference as deref

cdef class LocalAssembler:
    pass

cdef class RealIntegralOperatorLocalAssembler(LocalAssembler):
    pass

cdef class ComplexIntegralOperatorLocalAssembler(LocalAssembler):
    pass

cdef class LocalOperatorLocalAssembler(LocalAssembler):
    pass


cdef extern from "bempp/assembly/dense_global_block_assembler.hpp" namespace "Bempp":
    cdef shared_ptr[const c_DiscreteBoundaryOperator[RESULT]] c_assembleDenseBlock "Bempp::assembleDenseBlock"[BASIS,RESULT](
            int, int, int, int, const c_Space[BASIS]&, const c_Space[BASIS]&, const c_LocalAssemblerForIntegralOperators[RESULT]&, const c_ParameterList&) except +catch_exception

def assemble_dense_block_ext(rows, cols, Space domain not None, Space dual_to_range not None, assembler, ParameterList parameters not None):
    """Assemble a given subblock of a dense boundary operator."""

    cdef RealDiscreteBoundaryOperator real_discrete_operator = RealDiscreteBoundaryOperator()
    cdef ComplexDiscreteBoundaryOperator complex_discrete_operator = ComplexDiscreteBoundaryOperator()

    if rows == -1:
        row_start = 0
        row_end = dual_to_range.global_dof_count
    else:
        row_start = rows[0]
        row_end = rows[1]

    if cols == -1:
        col_start = 0
        col_end = dual_to_range.global_dof_count
    else:
        col_start = cols[0]
        col_end = cols[1]

    if isinstance(assembler, RealIntegralOperatorLocalAssembler):
        real_discrete_operator.impl_ = c_assembleDenseBlock[double,double](row_start, row_end, col_start, col_end,
                deref(dual_to_range.impl_), deref(domain.impl_), deref((<RealIntegralOperatorLocalAssembler> assembler).impl_), deref(parameters.impl_))
        return real_discrete_operator
    if isinstance(assembler, ComplexIntegralOperatorLocalAssembler):
        complex_discrete_operator.impl_ = c_assembleDenseBlock[double,complex_double](row_start, row_end, col_start, col_end,
                deref(dual_to_range.impl_), deref(domain.impl_), deref((<ComplexIntegralOperatorLocalAssembler> assembler).impl_), deref(parameters.impl_))
        return complex_discrete_operator
    raise ValueError("Unknown assembler type.")





