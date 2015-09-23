// Copyright (C) 2011 by the BEM++ Authors
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

#include "bempp/assembly/boundary_operator.hpp"
#include "bempp/assembly/evaluation_options.hpp"
#include "bempp/assembly/grid_function.hpp"
#include "bempp/assembly/l2_norm.hpp"
#include "bempp/assembly/numerical_quadrature_strategy.hpp"
#include "bempp/assembly/surface_normal_and_domain_index_dependent_function.hpp"
#include "bempp/assembly/discrete_boundary_operator.hpp"
#include "bempp/common/global_parameters.hpp"

#include "bempp/assembly/identity_operator.hpp"
#include "bempp/assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_single_layer_potential_operator.hpp"
#include "bempp/assembly/laplace_3d_double_layer_potential_operator.hpp"

#include "bempp/common/boost_make_shared_fwd.hpp"

#include "bempp/grid/grid.hpp"
#include "bempp/grid/grid_factory.hpp"

#include "bempp/space/piecewise_linear_continuous_scalar_space.hpp"
#include "bempp/space/piecewise_constant_scalar_space.hpp"


#include "bempp/common/eigen_support.hpp"

#include <iostream>
#include <fstream>

typedef double BFT; // basis function type
typedef double RT; // result type (type used to represent discrete operators)
typedef double CT; // coordinate type

class DirichletData
{
public:
    // Type of the function's values (e.g. float or std::complex<double>)
    typedef RT ValueType;
    // Type of coordinates (must be the "real part" of ValueType)
    typedef CT CoordinateType;

    // Number of components of the function's argument
    int argumentDimension() const { return 3; }
    // Number of components of the function's value
    int resultDimension() const { return 1; }

    // Evaluate the function at the point "point" and store result in
    // the array "result"
    inline void evaluate(const Eigen::Ref<Bempp::Vector<CoordinateType>>& point,
                         const Eigen::Ref<Bempp::Vector<ValueType>>& normal,
                         int domainIndex,
                         Eigen::Ref<Bempp::Vector<ValueType>> result) const {
        CoordinateType x = point(0), y = point(1), z = point(2);
        CoordinateType r = sqrt(point(0) * point(0) +
                point(1) * point(1) +
                point(2) * point(2));
        result(0) = 2 * x * z / (r * r * r * r * r) - y / (r * r * r);
    }
};

class ExactNeumannData
{
public:
    // Type of the function's values (e.g. float or std::complex<double>)
    typedef RT ValueType;
    // Type of coordinates (must be the "real part" of ValueType)
    typedef CT CoordinateType;

    // Number of components of the function's argument
    int argumentDimension() const { return 3; }
    // Number of components of the function's value
    int resultDimension() const { return 1; }

    // Evaluate the function at the point "point" and store result in
    // the array "result"
    inline void evaluate(const Eigen::Ref<Bempp::Vector<CoordinateType>>& point,
                         const Eigen::Ref<Bempp::Vector<CoordinateType>>& normal,
                         int domainIndex,
                         Eigen::Ref<Bempp::Vector<ValueType>> result) const {
        CoordinateType x = point(0), y = point(1), z = point(2);
        CoordinateType r = sqrt(point(0) * point(0) +
                point(1) * point(1) +
                point(2) * point(2));
        result(0) = -6 * x * z / (r * r * r * r * r * r) + 2 * y / (r * r * r * r);
    }
};

int main()
{
    // Import symbols from namespace Bempp to the global namespace

    using namespace Bempp;

    // Load mesh

    const char* meshFile = "/Users/betcke/development/bempp/meshes/sphere-h-0.05.msh";
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(params, meshFile);
    
    // Initialize the spaces

    PiecewiseLinearContinuousScalarSpace<BFT> pwiseLinears(grid);
    PiecewiseConstantScalarSpace<BFT> pwiseConstants(grid);

    // Define the Context object from default parameters
    
    auto parameters = GlobalParameters::parameterList();
    parameters.put("options.assembly.boundaryOperatorAssemblyType",std::string("hmat"));
    
    Context<BFT, RT> context(parameters);

    // Construct elementary operators

    BoundaryOperator<BFT, RT> slpOp =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseConstants),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseConstants));

    auto weak_form = slpOp.weakForm();
}
