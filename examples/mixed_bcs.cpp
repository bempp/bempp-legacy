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

// In this example we solve the Laplace equation with mixed boundary conditions
// (part Dirichlet, part Neumann).

#include "assembly/assembly_options.hpp"
#include "assembly/boundary_operator.hpp"
#include "assembly/blocked_boundary_operator.hpp"
#include "assembly/blocked_operator_structure.hpp"
#include "assembly/context.hpp"
#include "assembly/evaluation_options.hpp"
#include "assembly/grid_function.hpp"
#include "assembly/l2_norm.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "assembly/surface_normal_independent_function.hpp"

#include "assembly/identity_operator.hpp"
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_single_layer_potential_operator.hpp"
#include "assembly/laplace_3d_double_layer_potential_operator.hpp"

#include "common/boost_make_shared_fwd.hpp"

#include "grid/grid_view.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/entity.hpp"

#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"
#include "grid/grid_segment.hpp"

#include "linalg/default_iterative_solver.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"

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

class NeumannData
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

    const char* meshFile = "meshes/sphere-domains.msh";
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(params, meshFile);

    GridSegment segmentD = GridSegment::openDomain(*grid, 1);
    GridSegment segmentN = segmentD.complement();

    // Initialize the spaces

    // Despite variable names, only spaces of piecewise constants are used
    PiecewiseConstantScalarSpace<BFT> pwiseLinearsD(grid, segmentD);
    PiecewiseConstantScalarSpace<BFT> pwiseLinearsN(grid, segmentN);
    PiecewiseConstantScalarSpace<BFT> pwiseConstantsD(grid, segmentD);
    PiecewiseConstantScalarSpace<BFT> pwiseConstantsN(grid, segmentN);

    PiecewiseConstantScalarSpace<BFT> pwiseConstants(grid);
    PiecewiseConstantScalarSpace<BFT> pwiseLinears(grid);

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
    AcaOptions acaOptions; // Default parameters for ACA
    assemblyOptions.switchToAcaMode(acaOptions);

    // Create the assembly context

    Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

    // Construct elementary operators

    BoundaryOperator<BFT, RT> slpOpDD =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseConstantsD),
                make_shared_from_ref(pwiseLinearsD),
                make_shared_from_ref(pwiseConstantsD));
    BoundaryOperator<BFT, RT> slpOpDN =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseConstantsN),
                make_shared_from_ref(pwiseLinearsD),
                make_shared_from_ref(pwiseConstantsD));
    BoundaryOperator<BFT, RT> slpOpND =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseConstantsD),
                make_shared_from_ref(pwiseLinearsN),
                make_shared_from_ref(pwiseConstantsN));
    BoundaryOperator<BFT, RT> slpOpNN =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseConstantsN),
                make_shared_from_ref(pwiseLinearsN),
                make_shared_from_ref(pwiseConstantsN));

    BoundaryOperator<BFT, RT> dlpOpDD =
            laplace3dDoubleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinearsD),
                make_shared_from_ref(pwiseLinearsD),
                make_shared_from_ref(pwiseConstantsD));
    BoundaryOperator<BFT, RT> dlpOpDN =
            laplace3dDoubleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinearsN),
                make_shared_from_ref(pwiseLinearsD),
                make_shared_from_ref(pwiseConstantsD));
    BoundaryOperator<BFT, RT> dlpOpND =
            laplace3dDoubleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinearsD),
                make_shared_from_ref(pwiseLinearsN),
                make_shared_from_ref(pwiseConstantsN));
    BoundaryOperator<BFT, RT> dlpOpNN =
            laplace3dDoubleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinearsN),
                make_shared_from_ref(pwiseLinearsN),
                make_shared_from_ref(pwiseConstantsN));

    BoundaryOperator<BFT, RT> idOpDD =
            identityOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinearsD),
                make_shared_from_ref(pwiseLinearsD),
                make_shared_from_ref(pwiseConstantsD));
    BoundaryOperator<BFT, RT> idOpDN =
            identityOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinearsN),
                make_shared_from_ref(pwiseLinearsD),
                make_shared_from_ref(pwiseConstantsD));
    BoundaryOperator<BFT, RT> idOpND =
            identityOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinearsD),
                make_shared_from_ref(pwiseLinearsN),
                make_shared_from_ref(pwiseConstantsN));
    BoundaryOperator<BFT, RT> idOpNN =
            identityOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinearsN),
                make_shared_from_ref(pwiseLinearsN),
                make_shared_from_ref(pwiseConstantsN));

    // Build block operators

    BlockedOperatorStructure<BFT, RT> lftOpStruct;
    lftOpStruct.setBlock(0, 0, slpOpDD);
    lftOpStruct.setBlock(0, 1, 0.5 * idOpDN - dlpOpDN);
    lftOpStruct.setBlock(1, 0, slpOpND);
    lftOpStruct.setBlock(1, 1, 0.5 * idOpNN - dlpOpNN);
    BlockedBoundaryOperator<BFT, RT> lftOp(lftOpStruct);

    BlockedOperatorStructure<BFT, RT> rtOpStruct;
    rtOpStruct.setBlock(0, 0, -0.5 * idOpDD + dlpOpDD);
    rtOpStruct.setBlock(0, 1, -slpOpDN);
    rtOpStruct.setBlock(1, 0, -0.5 * idOpND + dlpOpND);
    rtOpStruct.setBlock(1, 1, -slpOpNN);
    BlockedBoundaryOperator<BFT, RT> rtOp(rtOpStruct);

    // Construct the grid functions representing the input data and the exact solutions

    GridFunction<BFT, RT> dirichletDataD(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinearsD),
                make_shared_from_ref(pwiseLinearsD),
                surfaceNormalIndependentFunction(DirichletData()));
    GridFunction<BFT, RT> exactDirichletDataN(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinearsN),
                make_shared_from_ref(pwiseLinearsN),
                surfaceNormalIndependentFunction(DirichletData()));

    GridFunction<BFT, RT> exactNeumannDataD(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseConstantsD),
                make_shared_from_ref(pwiseConstantsD),
                surfaceNormalIndependentFunction(NeumannData()));
    GridFunction<BFT, RT> neumannDataN(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseConstantsN),
                make_shared_from_ref(pwiseConstantsN),
                surfaceNormalIndependentFunction(NeumannData()));

    // Construct the right-hand-side grid function

    std::vector<GridFunction<BFT, RT> > rtArgs(2);
    rtArgs[0] = dirichletDataD;
    rtArgs[1] = neumannDataN;
    std::vector<GridFunction<BFT, RT> > rhs = rtOp * rtArgs;

    // Initialize the solver

    DefaultIterativeSolver<BFT, RT> solver(lftOp);
    solver.initializeSolver(defaultGmresParameterList(1e-8));

    // Solve the equation

    BlockedSolution<BFT, RT> solution = solver.solve(rhs);
    std::cout << solution.solverMessage() << std::endl;

    // Extract the solution in the form of grid functions
    // and export it in VTK format

    const GridFunction<BFT, RT>& neumannDataD = solution.gridFunction(0);
    const GridFunction<BFT, RT>& dirichletDataN = solution.gridFunction(1);

    dirichletDataD.exportToVtk(VtkWriter::CELL_DATA, "dD", "dD");
    dirichletDataN.exportToVtk(VtkWriter::CELL_DATA, "dN", "dN");
    neumannDataD.exportToVtk(VtkWriter::CELL_DATA, "nD", "nD");
    neumannDataN.exportToVtk(VtkWriter::CELL_DATA, "nN", "nN");
    exactDirichletDataN.exportToVtk(VtkWriter::CELL_DATA, "edN", "edN");
    exactNeumannDataD.exportToVtk(VtkWriter::CELL_DATA, "enD", "enD");

    // Construct "scattering operators"

    BoundaryOperator<BFT, RT> plD_to_pl =
            identityOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinearsD),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears));
    BoundaryOperator<BFT, RT> plN_to_pl =
            identityOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinearsN),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears));
    BoundaryOperator<BFT, RT> pcD_to_pc =
            identityOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseConstantsD),
                make_shared_from_ref(pwiseConstants),
                make_shared_from_ref(pwiseConstants));
    BoundaryOperator<BFT, RT> pcN_to_pc =
            identityOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseConstantsN),
                make_shared_from_ref(pwiseConstants),
                make_shared_from_ref(pwiseConstants));

    GridFunction<BFT, RT> dirichletData =
        plD_to_pl * dirichletDataD +
        plN_to_pl * exactDirichletDataN;
    GridFunction<BFT, RT> neumannData =
        pcD_to_pc * exactNeumannDataD +
        pcN_to_pc * neumannDataN;

    dirichletData.exportToVtk(VtkWriter::CELL_DATA, "dirichlet_data", "dirichlet_data");
    neumannData.exportToVtk(VtkWriter::CELL_DATA, "neumann_data", "neumann_data");

    CT absoluteError, relativeError;
    estimateL2Error(
                dirichletData, surfaceNormalIndependentFunction(DirichletData()),
                quadStrategy, absoluteError, relativeError);
    std::cout << "Relative L^2 error (Dirichlet): " << relativeError << std::endl;
    estimateL2Error(
                neumannData, surfaceNormalIndependentFunction(NeumannData()),
                quadStrategy, absoluteError, relativeError);
    std::cout << "Relative L^2 error (Neumann): " << relativeError << std::endl;
}
