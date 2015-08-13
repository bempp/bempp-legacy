from bempp_ext.utils cimport Matrix
from bempp_ext.utils cimport eigen_matrix_to_np_float64
from bempp_ext.utils cimport catch_exception
from bempp_ext.utils cimport shared_ptr, complex_double
from bempp_ext.utils cimport Connection, SlotInterface
from bempp_ext.grid.grid cimport Grid, c_Grid
from bempp_ext.grid.entity cimport Entity0, c_Entity
from bempp_ext.grid.codim_template cimport codim_zero
from libcpp cimport complex as ccomplex, bool as cbool
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "bempp/space/space.hpp":
    cdef cppclass c_Space "Bempp::Space" [BASIS]:
        c_Space(const shared_ptr[c_Grid]&)
        c_Space(const c_Space[BASIS]&)
        shared_ptr[const c_Grid] grid() const
        cbool spaceIsCompatible(const c_Space&)
        cbool is_same "isSame"(const c_Space&)
        int codomainDimension() const
        int domainDimension() const
        unsigned long globalDofCount() const
        unsigned long flatLocalDofCount() const 
        Connection connect(const SlotInterface&) const
        void getGlobalDofInterpolationPoints(Matrix[double]& points) const 
        void getNormalsAtGlobalDofInterpolationPoints(Matrix[double]& normals) const
        void getGlobalDofs(const c_Entity[codim_zero]&, vector[int]&, vector[double]&) const


cdef class Space:
    cdef:
        shared_ptr[const c_Space[double]] impl_

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

