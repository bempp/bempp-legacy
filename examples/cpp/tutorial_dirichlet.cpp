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

#include "bempp/assembly/assembly_options.hpp"
#include "bempp/assembly/boundary_operator.hpp"
#include "bempp/assembly/context.hpp"
#include "bempp/assembly/evaluation_options.hpp"
#include "bempp/assembly/grid_function.hpp"
#include "bempp/assembly/l2_norm.hpp"
#include "bempp/assembly/numerical_quadrature_strategy.hpp"
#include "bempp/assembly/surface_normal_independent_function.hpp"

#include "bempp/assembly/identity_operator.hpp"
#include "bempp/assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_single_layer_potential_operator.hpp"
#include "bempp/assembly/laplace_3d_double_layer_potential_operator.hpp"

#include "bempp/common/boost_make_shared_fwd.hpp"

#include "bempp/grid/grid.hpp"
#include "bempp/grid/grid_factory.hpp"

#include "bempp/linalg/default_iterative_solver.hpp"

#include "bempp/space/piecewise_linear_continuous_scalar_space.hpp"
#include "bempp/space/piecewise_constant_scalar_space.hpp"

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
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
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
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
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

    const char* meshFile = "../../../meshes/sphere-h-0.1.msh";
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(params, meshFile);
    

    // Initialize the spaces

    PiecewiseLinearContinuousScalarSpace<BFT> pwiseLinears(grid);
    PiecewiseConstantScalarSpace<BFT> pwiseConstants(grid);

    // Define the quadrature strategy

    AccuracyOptions accuracyOptions;
    // Increase by 2 the order of quadrature rule used to approximate
    // integrals of regular functions on pairs on elements
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    // Increase by 2 the order of quadrature rule used to approximate
    // integrals of regular functions on single elements
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);

    // Specify the assembly method. We want to use ACA

    AssemblyOptions assemblyOptions;
    // AcaOptions acaOptions; // Default parameters for ACA
    // assemblyOptions.switchToAcaMode(acaOptions);

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

    // Extract the solution in the form of a grid function
    // and export it in VTK format

    const GridFunction<BFT, RT>& solFun = solution.gridFunction();
    exportToVtk(solFun, VtkWriter::CELL_DATA, "Neumann_data", "solution");

    // Compare the numerical and analytical solution on the grid

    // GridFunction<BFT, RT> exactSolFun(
    //             make_shared_from_ref(context),
    //             make_shared_from_ref(pwiseConstants),
    //             make_shared_from_ref(pwiseConstants),
    //             surfaceNormalIndependentFunction(ExactNeumannData()));
    CT absoluteError, relativeError;
    estimateL2Error(
                solFun, surfaceNormalIndependentFunction(ExactNeumannData()),
                quadStrategy, absoluteError, relativeError);
    std::cout << "Relative L^2 error: " << relativeError << std::endl;

    // GridFunction<BFT, RT> diff = solFun - exactSolFun;
    // double relativeError = diff.L2Norm() / exactSolFun.L2Norm();
    // std::cout << "Relative L^2 error: " << relativeError << std::endl;

    // Prepare to evaluate the solution on an annulus outside the sphere

    // Create potential operators

    Laplace3dSingleLayerPotentialOperator<BFT, RT> slPotOp;
    Laplace3dDoubleLayerPotentialOperator<BFT, RT> dlPotOp;

    // Construct the array 'evaluationPoints' containing the coordinates
    // of points where the solution should be evaluated

    const int rCount = 51;
    const int thetaCount = 361;
    const CT minTheta = 0., maxTheta = 2. * M_PI;
    const CT minR = 1., maxR = 2.;
    const int dimWorld = 3;
    arma::Mat<CT> evaluationPoints(dimWorld, rCount * thetaCount);
    for (int iTheta = 0; iTheta < thetaCount; ++iTheta) {
        CT theta = minTheta + (maxTheta - minTheta) *
            iTheta / (thetaCount - 1);
        for (int iR = 0; iR < rCount; ++iR) {
            CT r = minR + (maxR - minR) * iR / (rCount - 1);
            evaluationPoints(0, iR + iTheta * rCount) = r * cos(theta); // x
            evaluationPoints(1, iR + iTheta * rCount) = r * sin(theta); // y
            evaluationPoints(2, iR + iTheta * rCount) = 0.;             // z
        }
    }

    // Use the Green's representation formula to evaluate the solution

    EvaluationOptions evaluationOptions;

    arma::Mat<RT> field =
        -slPotOp.evaluateAtPoints(solFun, evaluationPoints,
                                  quadStrategy, evaluationOptions) +
         dlPotOp.evaluateAtPoints(dirichletData, evaluationPoints,
                                  quadStrategy, evaluationOptions);

    // Export the solution into text file

    std::ofstream out("solution.txt");
    out << "# x y z u\n";
    for (int i = 0; i < rCount * thetaCount; ++i)
        out << evaluationPoints(0, i) << ' '
            << evaluationPoints(1, i) << ' '
            << evaluationPoints(2, i) << ' '
            << field(0, i) << '\n';
}
