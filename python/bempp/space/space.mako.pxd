<%
from data_types import dtypes, scalar_cython_type, real_cython_type

from space import spaces

def declare_class(text):
    return 'c_{0} "Bempp::{0}"[BASIS]'.format(text)
%>
from bempp.utils cimport Matrix
from bempp.utils.eigen cimport eigen_matrix_to_np_float64
from bempp.utils cimport catch_exception
from libcpp cimport complex as ccomplex, bool as cbool
from libcpp.string cimport string
from libcpp.vector cimport vector
from bempp.utils cimport shared_ptr, complex_float,complex_double
from bempp.utils.signal_slot_interface cimport Connection, SlotInterface
from bempp.grid.grid cimport Grid, c_Grid
from bempp.grid.entity cimport Entity0, c_Entity
from bempp.grid.codim_template cimport codim_zero

cdef extern from "bempp/space/py_space_variants.hpp":
    cdef cppclass c_Space "Bempp::Space" [BASIS]:
        c_Space(const shared_ptr[c_Grid]&)
        c_Space(const c_Space[BASIS]&)
        shared_ptr[const c_Grid] grid() const

        
# Declares complex type explicitly.
# Cython 0.20 will fail if templates are nested more than three-deep,
# as in shared_ptr[ c_Space[ complex[float] ] ]
cdef extern from "bempp/utils/py_types.hpp":
% for ctype in dtypes.values():
%     if 'complex'  in ctype:
    ctypedef struct ${ctype}
%     endif
% endfor



% for class_name, description in spaces.items():
cdef extern from "${description['header']}":
    cdef cppclass ${class_name | declare_class}:
%   if description['implementation'] == 'grid_only':
        ${'c_' + class_name}(const shared_ptr[c_Grid] &_grid)
%   elif description['implementation'] == 'polynomial':
        ${'c_' + class_name}(const shared_ptr[c_Grid] &_grid, int order)
%   endif
% endfor

cdef extern from "bempp/space/py_space_variants.hpp" namespace "Bempp":
    cdef cppclass SpaceVariants:
        SpaceVariants()
        void set[T](const shared_ptr[T] &_in)
        void set[T](const shared_ptr[const T] &_in)
        # void reset[T](const shared_ptr[T] &_in)
        string dtype() const
        shared_ptr[const c_Grid] grid() const
        cbool isCompatible(const SpaceVariants&)
        cbool is_same "isSame"(const SpaceVariants&)
        int codomainDimension() const
        int domainDimension() const
        unsigned long globalDofCount() const
        unsigned long flatLocalDofCount() const 
        Connection connect(const SlotInterface&) const

    cdef shared_ptr[c_Space[BASIS]] _py_get_space_ptr[BASIS](const SpaceVariants& space_variant)
% for pybasis,cybasis in dtypes.items():
    cdef Matrix[${real_cython_type(cybasis)}] _py_space_get_global_dof_interp_points_${pybasis} "Bempp::_py_space_get_global_dof_interp_points<${cybasis}>"(const SpaceVariants& space_variant) except +catch_exception
    cdef Matrix[${real_cython_type(cybasis)}] _py_space_get_global_dof_normals_${pybasis} "Bempp::_py_space_get_global_dof_normals<${cybasis}>"(const SpaceVariants& space_variant) except +catch_exception
    cdef void _py_space_get_global_dofs_${pybasis} "Bempp::_py_space_get_global_dofs<${cybasis}>"(const SpaceVariants& space_variant, const c_Entity[codim_zero]&, vector[int]&, vector[${cybasis}]&) except +catch_exception
% endfor

cdef class Space:
    cdef:
        SpaceVariants impl_
        unsigned int _order
    cpdef cbool is_compatible(self, Space other)

# Define all possible spaces

cdef extern from "bempp/space/piecewise_constant_scalar_space.hpp" namespace "Bempp":
    cdef shared_ptr[c_Space[T]] adaptivePiecewiseConstantScalarSpace[T](const shared_ptr[c_Grid]& grid)
    cdef shared_ptr[c_Space[T]] adaptivePiecewiseConstantScalarSpace[T](const shared_ptr[c_Grid]& grid, vector[int] domains, cbool closed)
cdef extern from "bempp/space/piecewise_linear_continuous_scalar_space.hpp" namespace "Bempp":
    cdef shared_ptr[c_Space[T]] adaptivePiecewiseLinearContinuousScalarSpace[T](const shared_ptr[c_Grid]& grid)
    cdef shared_ptr[c_Space[T]] adaptivePiecewiseLinearContinuousScalarSpace[T](const shared_ptr[c_Grid]& grid, vector[int] domains, cbool closed)
cdef extern from "bempp/space/piecewise_linear_discontinuous_scalar_space.hpp" namespace "Bempp":
    cdef shared_ptr[c_Space[T]] adaptivePiecewiseLinearDiscontinuousScalarSpace[T](const shared_ptr[c_Grid]& grid)
    cdef shared_ptr[c_Space[T]] adaptivePiecewiseLinearDiscontinuousScalarSpace[T](const shared_ptr[c_Grid]& grid, vector[int] domains, cbool closed)
cdef extern from "bempp/space/piecewise_polynomial_continuous_scalar_space.hpp" namespace "Bempp":
    cdef shared_ptr[c_Space[T]] adaptivePiecewisePolynomialContinuousScalarSpace[T](const shared_ptr[c_Grid]& grid, int order)
    cdef shared_ptr[c_Space[T]] adaptivePiecewisePolynomialContinuousScalarSpace[T](const shared_ptr[c_Grid]& grid, int order, vector[int] domains, cbool closed)
cdef extern from "bempp/space/piecewise_polynomial_discontinuous_scalar_space.hpp" namespace "Bempp":
    cdef shared_ptr[c_Space[T]] adaptivePiecewisePolynomialDiscontinuousScalarSpace[T](const shared_ptr[c_Grid]& grid, int order)
    cdef shared_ptr[c_Space[T]] adaptivePiecewisePolynomialDiscontinuousScalarSpace[T](const shared_ptr[c_Grid]& grid, int order, vector[int] domains, cbool closed)
cdef extern from "bempp/space/raviart_thomas_0_vector_space.hpp" namespace "Bempp":
    cdef shared_ptr[c_Space[T]] adaptiveRaviartThomas0VectorSpace[T](const shared_ptr[c_Grid]& grid)
    cdef shared_ptr[c_Space[T]] adaptiveRaviartThomas0VectorSpace[T](const shared_ptr[c_Grid]& grid, vector[int] domains, cbool closed)

