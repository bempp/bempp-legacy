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

#include "assembly/assembly_options.hpp"
#include "assembly/boundary_operator.hpp"
#include "assembly/context.hpp"
#include "assembly/grid_function.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "assembly/surface_normal_independent_function.hpp"

#include "assembly/identity_operator.hpp"
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_double_layer_boundary_operator.hpp"

#include "common/boost_make_shared_fwd.hpp"

#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"

#include "linalg/default_iterative_solver.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"

#include <iostream>

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
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        CoordinateType x = point(0), y = point(1), z = point(2);
        CoordinateType r = arma::norm(point, 2);
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
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        CoordinateType x = point(0), y = point(1), z = point(2);
        CoordinateType r = arma::norm(point, 2);
        result(0) = -6 * x * z / (r * r * r * r * r * r) + 2 * y / (r * r * r * r);
    }
};

int main()
{
    // Import symbols from namespace Bempp to the global namespace

    using namespace Bempp;

    // Load mesh

    const char* meshFile = "meshes/sphere-644.msh";
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    std::auto_ptr<Grid> grid = GridFactory::importGmshGrid(params, meshFile);

    // Initialize the spaces

    PiecewiseLinearContinuousScalarSpace<BFT> pwiseLinears(*grid);
    PiecewiseConstantScalarSpace<BFT> pwiseConstants(*grid);

    pwiseLinears.assignDofs();
    pwiseConstants.assignDofs();

    // Define the quadrature strategy

    NumericalQuadratureStrategy<BFT, RT> quadStrategy;

    // Specify the assembly method. We want to use ACA

    AssemblyOptions assemblyOptions;
    AcaOptions acaOptions; // Default parameters for ACA
    assemblyOptions.switchToAcaMode(acaOptions);

    // Create the assembly context

    Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

    // Construct elementary operators

    BoundaryOperator<BFT, RT> slpOp =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseConstants),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseConstants));
    BoundaryOperator<BFT, RT> dlpOp =
            laplace3dDoubleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseConstants));
    BoundaryOperator<BFT, RT> idOp =
            identityOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseConstants));

    // Form the right-hand side sum

    BoundaryOperator<BFT, RT> rhsOp = -0.5 * idOp + dlpOp;

    // Construct the grid function representing the (input) Dirichlet data

    GridFunction<BFT, RT> dirichletData(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                surfaceNormalIndependentFunction(DirichletData()));

    // Construct the right-hand-side grid function 

    GridFunction<BFT, RT> rhs = rhsOp * dirichletData;

    // Initialize the solver

    DefaultIterativeSolver<BFT, RT> solver(slpOp);
    solver.initializeSolver(defaultGmresParameterList(1e-5));

    // Solve the equation

    Solution<BFT, RT> solution = solver.solve(rhs);
    std::cout << solution.solverMessage() << std::endl;

    // Extract the solution in the form of a grid function and export it in VTK format

    const GridFunction<BFT, RT>& solFun = solution.gridFunction();
    solFun.exportToVtk(VtkWriter::CELL_DATA, "Neumann_data", "solution");

    // Compare the numerical and analytical solution

    GridFunction<BFT, RT> exactSolFun(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseConstants),
                make_shared_from_ref(pwiseConstants),
                surfaceNormalIndependentFunction(ExactNeumannData()));
    GridFunction<BFT, RT> diff = solFun - exactSolFun;
    double relativeError = diff.L2Norm() / exactSolFun.L2Norm();
    std::cout << "Relative L^2 error: " << relativeError << std::endl;
}
