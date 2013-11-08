// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// Construct a Swig module
%module(directors="1") core
%{
#define SWIG_FILE_WITH_INIT

#include <dune/common/exceptions.hh>
#include <complex>
%}

%include "numpy.i"

%include "std_vector.i"
%include "std_string.i"

namespace std {

    %template(IntVector) vector<int>;
    %template(DoubleVector) vector<double>;
    %template(VectorDoubleVector) vector<vector<double> >;
    %template(StringVector) vector<std::string>;

}

%include "config.i"

// Import docstring macros
%include "docstrings.i"

// Define commonly used typemaps
%include "armadillo.i"
%include "auto_ptr.i"
%include "exception.i"
%include "bempp_shared_ptr.i"
%include "std_string.i"
%include "std_complex.i"

// Useful macros
%include "macros.i"

// Some useful Python routines
%include "py_defs.i"

// Import macros for explicit template instantiations
%include "template_instantiations_basis.i"
%include "template_instantiations_basis_result.i"
%include "template_instantiations_basis_result_geometry_factory.i"
%include "template_instantiations_basis_kernel_result.i"
%include "template_instantiations_value.i"

// Setup a handler for C++ exceptions
%exception {
    try {
        $action
    }
    catch (const Dune::Exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what().c_str());
    }
    catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}

// Declare all necessary auto_ptr typemaps

// Grid
AUTO_PTR_TYPEMAPS(Bempp::Grid)
AUTO_PTR_TYPEMAPS(Bempp::GridView)
AUTO_PTR_TYPEMAPS(Bempp::EntityPointer<0>)
AUTO_PTR_TYPEMAPS(Bempp::EntityPointer<1>)
AUTO_PTR_TYPEMAPS(Bempp::EntityPointer<2>)
AUTO_PTR_TYPEMAPS(Bempp::EntityPointer<3>)
AUTO_PTR_TYPEMAPS(Bempp::EntityIterator<0>)
AUTO_PTR_TYPEMAPS(Bempp::EntityIterator<1>)
AUTO_PTR_TYPEMAPS(Bempp::EntityIterator<2>)
AUTO_PTR_TYPEMAPS(Bempp::EntityIterator<3>)
AUTO_PTR_TYPEMAPS(Bempp::Geometry)
AUTO_PTR_TYPEMAPS(Bempp::VtkWriter)

AUTO_PTR_TYPEMAPS_FOR_CLASS_TEMPLATED_ON_RESULT(Bempp::DiscreteBoundaryOperator)
AUTO_PTR_TYPEMAPS_FOR_CLASS_TEMPLATED_ON_RESULT(Bempp::InterpolatedFunction)

// End of auto_ptr typemaps

%{
#include "common/deprecated.hpp"
%}
#define BEMPP_DEPRECATED

%feature("autodoc", 0);

// Make commonly used typedefs known to Swig
%inline %{
    namespace Bempp
    {
        typedef double ctype;
    }
%}
%include "grid/geometry_type.hpp"

// Wrap Bempp components

// Common
%include "common/scalar_traits.i"
// Grid
%include "grid/grid_parameters.i"
%include "grid/geometry.i"
%include "grid/geometry_type.i"
%include "grid/entity.i"
%include "grid/entity_pointer.i"
%include "grid/entity_iterator.i"
%include "grid/grid.i"
%include "grid/grid_segment.i"
%include "grid/id_set.i"
%include "grid/grid_view.i"
%include "grid/index_set.i"
%include "grid/vtk_writer.i"
%include "grid/grid_factory.i"
%include "grid/geometry_factory.i"

%feature("autodoc", 2);

// Fiber
%include "fiber/opencl_options.i"
%include "fiber/parallelization_options.i"
%include "fiber/quadrature_options.i"
%include "fiber/accuracy_options.i"
%include "fiber/quadrature_strategy.i"
%include "fiber/verbosity_level.i"

// Space
%include "space/space.i"
%include "space/piecewise_constant_scalar_space.i"
%include "space/piecewise_constant_scalar_space_barycentric.i"
%include "space/piecewise_constant_dual_grid_scalar_space.i"
%include "space/piecewise_linear_continuous_scalar_space.i"
%include "space/piecewise_linear_discontinuous_scalar_space.i"
%include "space/piecewise_polynomial_continuous_scalar_space.i"
%include "space/piecewise_polynomial_discontinuous_scalar_space.i"
%include "space/piecewise_linear_continuous_scalar_space_barycentric.i"
%include "space/piecewise_linear_discontinuous_scalar_space_barycentric.i"
%include "space/raviart_thomas_0_vector_space.i"
%include "space/unit_scalar_space.i"

// Assembly
%include "assembly/aca_options.i"
%include "assembly/assembly_options.i"
%include "assembly/numerical_quadrature_strategy.i"
%include "assembly/transposition_mode.i"
%include "assembly/symmetry.i"
%include "assembly/python_domain_index_dependent_functor.i"
%include "assembly/python_surface_normal_and_domain_index_dependent_functor.i"
%include "assembly/python_surface_normal_dependent_functor.i"
%include "assembly/python_surface_normal_independent_functor.i"
%include "assembly/discrete_boundary_operator.i"
%include "assembly/context.i"
%include "assembly/grid_function.i"
%include "assembly/evaluation_options.i"
%include "assembly/l2_norm.i"
%include "assembly/abstract_boundary_operator.i"
%include "assembly/boundary_operator.i"
%include "assembly/laplace_3d_operators.i"
%include "assembly/helmholtz_3d_operators.i"
%include "assembly/modified_helmholtz_3d_operators.i"
%include "assembly/maxwell_3d_operators.i"
%include "assembly/identity_operator.i"
%include "assembly/null_operator.i"
%include "assembly/laplace_beltrami_3d_operator.i"
%include "assembly/potential_operator.i"
%include "assembly/assembled_potential_operator.i"
%include "assembly/laplace_3d_potential_operators.i"
%include "assembly/helmholtz_3d_potential_operators.i"
%include "assembly/modified_helmholtz_3d_potential_operators.i"
%include "assembly/maxwell_3d_potential_operators.i"
%include "assembly/blocked_operator_structure.i"
%include "assembly/blocked_boundary_operator.i"
%include "assembly/discrete_aca_boundary_operator.i"
%include "assembly/discrete_dense_boundary_operator.i"
%include "assembly/discrete_inverse_sparse_boundary_operator.i"

// Linear algebra
%include "linalg/parameter_list.i"
%include "linalg/solution_base.i"
%include "linalg/solution.i"
%include "linalg/blocked_solution.i"
%include "linalg/solver.i"
%include "linalg/preconditioner.i"
%include "linalg/default_iterative_solver.i"
%include "linalg/default_direct_solver.i"

// IO

%include "io/gmsh.i"
