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

// This script solves the Maxwell equations in the region exterior to a bounded
// object, with Dirichlet boundary conditions given by the exact solution
// (satisfying the Silver-Mueller radiation conditions)
//
//     \vec u(\vec x) = h_1^{(1)}(k r) sin(theta) \hat phi,
//
// where (r, theta, phi) are the radial, zenith angle and azimuthal spherical
// coordinates, h_1^{(1)}(r) is the spherical Hankel function of the first kind
// and order 1 and \hat phi is the unit vector oriented along d(\vec x)/d\phi.
//
// Many options can be controlled from the command line (currently the processing
// of command-line options is ugly -- just positional arguments; in future,
// proper named arguments will be added). Example invocation:
//
//     ./maxwell_dirichlet sphere-ico-3.msh -1 1e-5 1e-8 0 2
//
// (see the usage message in the source code for the description of specific
// arguments).

#include "assembly/aca_approximate_lu_inverse.hpp"
#include "assembly/assembly_options.hpp"
#include "assembly/boundary_operator.hpp"
#include "assembly/context.hpp"
#include "assembly/discrete_aca_boundary_operator.hpp"
#include "assembly/evaluation_options.hpp"
#include "assembly/grid_function.hpp"
#include "assembly/l2_norm.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "assembly/surface_normal_dependent_function.hpp"
#include "assembly/surface_normal_independent_function.hpp"

#include "assembly/maxwell_3d_identity_operator.hpp"
#include "assembly/maxwell_3d_single_layer_boundary_operator.hpp"
#include "assembly/maxwell_3d_double_layer_boundary_operator.hpp"
#include "assembly/maxwell_3d_single_layer_potential_operator.hpp"
#include "assembly/maxwell_3d_double_layer_potential_operator.hpp"

#include "common/armadillo_fwd.hpp"
#include "common/boost_make_shared_fwd.hpp"
#include "common/scalar_traits.hpp"

#include "grid/grid_factory.hpp"
#include "grid/grid.hpp"

#include "linalg/preconditioner.hpp"
#include "linalg/default_iterative_solver.hpp"
#include "linalg/default_direct_solver.hpp"

#include "space/piecewise_linear_normally_continuous_vector_space.hpp"

#include <cmath>
#include <iostream>
#include <memory>

using namespace Bempp;

typedef double CT; // coordinate type
typedef CT BFT; // basis function type
typedef std::complex<CT> RT; // result type (type used to represent elements of
                             // discrete operators)

const CT k = 2; // wave number

const CT source = 0.1; // the source will be located at the point
                       // (source, source, source)

// Return the cross product of the vectors a and b
template <typename T, typename U>
arma::Col<typename Fiber::Coercion<T, U>::Type> crossProduct(
        const arma::Col<T>& a, const arma::Col<U>& b)
{
    assert(a.n_rows == 3 && b.n_rows == 3);
    arma::Col<typename Fiber::Coercion<T, U>::Type> result(3);
    result(0) = a(1) * b(2) - a(2) * b(1);
    result(1) = a(2) * b(0) - a(0) * b(2);
    result(2) = a(0) * b(1) - a(1) * b(0);
    return result;
}

class DirichletData
{
public:
    // Type of the function's values (e.g. float or std::complex<double>)
    typedef RT ValueType;
    // Type of coordinates (must be the "real part" of ValueType)
    typedef CT CoordinateType;

    // Number of components of the function's argument
    int argumentDimension() const { return 3; }
    // Number of components of the function's result
    int resultDimension() const { return 3; }

    // Evaluate the function at the point "point" and store result in
    // the array "result"
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         const arma::Col<CoordinateType>& normal,
                         arma::Col<ValueType>& result) const {
        arma::Col<RT> field(3);
        CT x = point(0) - source, y = point(1) - source, z = point(2) - source;
        CT r = sqrt(x * x + y * y + z * z);
        const RT I(0., 1.);
        CT kr = k * r;
        RT h1kr = (-I - kr) * exp(I * kr) / (kr * kr);
        RT scale = h1kr / r;
        field(0) = -y * scale;
        field(1) = x * scale;
        field(2) = 0.;
        result = crossProduct(field, normal);
    }
};

class ExactNeumannData
{
public:
    typedef RT ValueType;
    typedef CT CoordinateType;

    int argumentDimension() const { return 3; }
    int resultDimension() const { return 3; }

    inline void evaluate(const arma::Col<CoordinateType>& point,
                         const arma::Col<CoordinateType>& normal,
                         arma::Col<ValueType>& result) const {
        arma::Col<RT> curl(3);
        CT x = point(0) - source, y = point(1) - source, z = point(2) - source;
        CT r = sqrt(x * x + y * y + z * z);
        const RT I(0., 1.);
        CT kr = k * r;
        RT h1kr = (-I - kr) * exp(I * kr) / (kr * kr);
        RT h1kr_deriv =
                (1. + I - I * kr) * (1. + I + kr) * exp(I * kr) / (kr * kr * r);
        RT xy_factor = (h1kr - r * h1kr_deriv) / (r * r * r);
        curl(0) = x * z * xy_factor;
        curl(1) = y * z * xy_factor;
        curl(2) = ((x*x + y*y + 2*z*z) * h1kr + r * (x*x + y*y) * h1kr_deriv) /
                (r * r * r);
        result = crossProduct(curl, normal) / (I * k);
    }
};

class ExactSolution
{
public:
    typedef RT ValueType;
    typedef CT CoordinateType;

    int argumentDimension() const { return 3; }
    int resultDimension() const { return 3; }

    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        arma::Col<RT> field(3);
        CT x = point(0) - source, y = point(1) - source, z = point(2) - source;
        CT r = sqrt(x * x + y * y + z * z);
        const RT I(0., 1.);
        CT kr = k * r;
        RT h1kr = (-I - k * r) * exp(I * k * r) / (kr * kr);
        RT scale = h1kr / r;
        field(0) = -y * scale;
        field(1) = x * scale;
        field(2) = 0.;
        result = field;
    }
};

int main(int argc, char* argv[])
{
    // Process command-line args

    if (argc < 7 || argc % 2 != 1) {
        std::cout << "Solve a Maxwell Dirichlet problem in an exterior domain.\n"
                  << "Usage: " << argv[0]
                  << " <mesh_file> <n_threads> <aca_eps> <solver_tol>"
                  << " <singular_order_increment>"
                  << " [<regular_order_increment_1> <min_relative_distance_1>]"
                  << " [<regular_order_increment_2> <min_relative_distance_2>]"
                  << " [...] <regular_order_increment_n>"
                  << std::endl;
        return 1;
    }
    int maxThreadCount = atoi(argv[2]);
    double acaEps = atof(argv[3]);
    double solverTol = atof(argv[4]);
    int singOrderIncrement = atoi(argv[5]);
    if (acaEps > 1. || acaEps < 0.) {
        std::cout << "Invalid aca_eps: " << acaEps << std::endl;
        return 1;
    }
    if (solverTol > 1. || solverTol < 0.) {
        std::cout << "Invalid solver_tol: " << solverTol << std::endl;
        return 1;
    }

    AccuracyOptionsEx accuracyOptions;
    std::vector<double> maxNormDists;
    std::vector<int> orderIncrements;
    for (int i = 6; i < argc - 1; i += 2) {
        orderIncrements.push_back(atoi(argv[i]));
        maxNormDists.push_back(atof(argv[i + 1]));
    }
    orderIncrements.push_back(atoi(argv[argc - 1]));
    accuracyOptions.setDoubleRegular(maxNormDists, orderIncrements);
    accuracyOptions.setDoubleSingular(singOrderIncrement);
    accuracyOptions.setSingleRegular(2);

    // Load mesh

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(params, argv[1]);

    // Initialize the space

    PiecewiseLinearNormallyContinuousVectorSpace<BFT> HdivSpace(grid);

    // Set assembly mode and options

    AssemblyOptions assemblyOptions;
    assemblyOptions.setMaxThreadCount(maxThreadCount);
    if (acaEps > 0.) {
        AcaOptions acaOptions;
        acaOptions.eps = acaEps;
        assemblyOptions.switchToAcaMode(acaOptions);
    }

    NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);
    Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

    // Construct operators

    BoundaryOperator<BFT, RT> slpOp = maxwell3dSingleLayerBoundaryOperator<BFT>(
                make_shared_from_ref(context),
                make_shared_from_ref(HdivSpace),
                make_shared_from_ref(HdivSpace),
                make_shared_from_ref(HdivSpace),
                k,
                "SLP");
    BoundaryOperator<BFT, RT> dlpOp = maxwell3dDoubleLayerBoundaryOperator<BFT>(
                make_shared_from_ref(context),
                make_shared_from_ref(HdivSpace),
                make_shared_from_ref(HdivSpace),
                make_shared_from_ref(HdivSpace),
                k,
                "DLP");
    BoundaryOperator<BFT, RT> idOp = maxwell3dIdentityOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(HdivSpace),
                make_shared_from_ref(HdivSpace),
                make_shared_from_ref(HdivSpace),
                "Id");

    // Construct a grid function representing the Dirichlet data

    GridFunction<BFT, RT> dirichletData(
                make_shared_from_ref(context),
                make_shared_from_ref(HdivSpace),
                make_shared_from_ref(HdivSpace),
                surfaceNormalDependentFunction(DirichletData()));

    dirichletData.exportToVtk(VtkWriter::CELL_DATA, "Dirichlet_data",
                              "input_dirichlet_data_cell");
    dirichletData.exportToVtk(VtkWriter::VERTEX_DATA, "Dirichlet_data",
                              "input_dirichlet_data_vertex");

    // Construct a grid function representing the right-hand side

    GridFunction<BFT, RT> rhs = -((0.5 * idOp + dlpOp) * dirichletData);

    // Solve the equation

    Solution<BFT, RT> solution;
    if (solverTol > 0.) {
        DefaultIterativeSolver<BFT, RT> solver(slpOp);

        AcaApproximateLuInverse<RT> slpOpLu(
                     DiscreteAcaBoundaryOperator<RT>::castToAca(
                         *slpOp.weakForm()),
                         /* LU factorisation accuracy */ 0.01);
        Preconditioner<RT> prec =
                     discreteOperatorToPreconditioner<RT>(
                         make_shared_from_ref(slpOpLu));

        solver.initializeSolver(defaultGmresParameterList(solverTol, 10000), prec);
        solution = solver.solve(rhs);
    } else {
        DefaultDirectSolver<BFT, RT> solver(slpOp);
        solution = solver.solve(rhs);
    }
    std::cout << solution.solverMessage() << std::endl;

    // Extract the solution

    const GridFunction<BFT, RT>& neumannData = solution.gridFunction();

    neumannData.exportToVtk(VtkWriter::CELL_DATA, "Neumann_data",
                        "calculated_neumann_data_cell");
    neumannData.exportToVtk(VtkWriter::VERTEX_DATA, "Neumann_data",
                        "calculated_neumann_data_vertex");

    // Compare the solution against the analytical result
    GridFunction<BFT, RT> exactNeumannData(
                make_shared_from_ref(context),
                make_shared_from_ref(HdivSpace),
                make_shared_from_ref(HdivSpace),
                surfaceNormalDependentFunction(ExactNeumannData()));
    exactNeumannData.exportToVtk(VtkWriter::CELL_DATA, "Neumann_data",
                                 "exact_neumann_data_cell");
    exactNeumannData.exportToVtk(VtkWriter::VERTEX_DATA, "Neumann_data",
                                 "exact_neumann_data_vertex");
    EvaluationOptions evaluationOptions;
    CT absoluteError = L2NormOfDifference(
                neumannData, surfaceNormalDependentFunction(ExactNeumannData()),
                quadStrategy, evaluationOptions);
    CT exactSolNorm = L2NormOfDifference(
                0. * neumannData, surfaceNormalDependentFunction(ExactNeumannData()),
                quadStrategy, evaluationOptions);
    std::cout << "Relative L^2 error: " << absoluteError / exactSolNorm
              << "\nAbsolute L^2 error: " << absoluteError << std::endl;

    // Evaluate the solution at a few points

    Maxwell3dSingleLayerPotentialOperator<BFT> slpPotOp(k);
    Maxwell3dDoubleLayerPotentialOperator<BFT> dlpPotOp(k);

    const int dimWorld = 3;
    const int pointCount = 3;
    arma::Mat<CT> points(dimWorld, pointCount * pointCount);
    for (int i = 0; i < pointCount; ++i)
        for (int j = 0; j < pointCount; ++j) {
            points(0, i * pointCount + j) = 3. + i / CT(pointCount - 1);
            points(1, i * pointCount + j) = 2. + j / CT(pointCount - 1);
            points(2, i * pointCount + j) = 0.5;
        }

    arma::Mat<RT> slpContrib = slpPotOp.evaluateAtPoints(
        neumannData, points, quadStrategy, evaluationOptions);
    arma::Mat<RT> dlpContrib = dlpPotOp.evaluateAtPoints(
        dirichletData, points, quadStrategy, evaluationOptions);
    arma::Mat<RT> values = -slpContrib - dlpContrib;

    // Evaluate the analytical solution at the same points

    ExactSolution exactSolution;
    arma::Mat<RT> exactValues(dimWorld, pointCount * pointCount);
    for (int i = 0; i < pointCount * pointCount; ++i) {
        arma::Col<RT> activeResultColumn = exactValues.unsafe_col(i);
        exactSolution.evaluate(points.unsafe_col(i), activeResultColumn);
    }

    // Compare the numerical and analytical solutions

    std::cout << "Numerical | analytical\n";
    for (int i = 0; i < pointCount * pointCount; ++i)
        std::cout << values(0, i) << " "
                  << values(1, i) << " "
                  << values(2, i) << " | "
                  << exactValues(0, i) << " "
                  << exactValues(1, i) << " "
                  << exactValues(2, i) << "\n";
}
