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

#include <armadillo>
#include <cmath>
#include <iostream>
#include <memory>

#include "meshes.hpp"

#include "assembly/assembly_options.hpp"
#include "assembly/discrete_scalar_valued_linear_operator.hpp"
#include "assembly/evaluation_options.hpp"
#include "assembly/grid_function.hpp"
#include "assembly/interpolated_function.hpp"

#include "assembly/identity_operator.hpp"
#include "assembly/single_layer_potential_3d.hpp"
#include "assembly/double_layer_potential_3d.hpp"
#include "assembly/adjoint_double_layer_potential_3d.hpp"
#include "assembly/hypersingular_operator_3d.hpp"

#include "fiber/standard_local_assembler_factory_for_operators_on_surfaces.hpp"

#include "grid/geometry.hpp"
#include "grid/geometry_factory.hpp"
#include "grid/grid_factory.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"

using namespace Bempp;

/**
  This script solves the Laplace equation

    \nabla^2 V = 0

  outside a sphere with radius 1, with the Dirichlet boundary condition

    V(x) = 1

  imposed on the surface of the sphere and the condition

    lim_{|x|\to\infty} V(x) = 0

  imposed at infinity. The analytical solution is

    V(x) = 1/|x|,

  where |x| is the distance from the centre of the sphere. Thus, the derivative
  of V on the surface in the (outer) normal direction is -1.
  */

int main()
{
    const MeshVariant meshVariant = SPHERICAL_SHELL_INNER_SURFACE;

    std::auto_ptr<Grid> grid = loadMesh(meshVariant);
    dumpElementList(grid.get());

    PiecewiseLinearContinuousScalarSpace<double> space(*grid);
    space.assignDofs();

    AssemblyOptions assemblyOptions;
    assemblyOptions.switchToDense();

    assemblyOptions.switchToTbb();
    assemblyOptions.setSingularIntegralCaching(AssemblyOptions::YES);

    EvaluationOptions evaluationOptions;
    evaluationOptions.switchToTbb();

    Fiber::OpenClOptions openClOptions;
    //    assemblyOptions.switchToOpenCl(openClOptions);

    Fiber::AccuracyOptions accuracyOptions; // default
    Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory>
            factory(accuracyOptions);

    typedef std::auto_ptr<LinearOperator<double> > LinearOperatorPtr;

    arma::Col<double> solutionCoefficients;

    SingleLayerPotential3D<double> slp;
    DoubleLayerPotential3D<double> dlp;
    IdentityOperator<double> id;

    typedef DiscreteScalarValuedLinearOperator<double> DiscreteLinearOperator;
    typedef std::auto_ptr<DiscreteLinearOperator> DiscreteLinearOperatorPtr;
    DiscreteLinearOperatorPtr discreteSlp =
            slp.assembleWeakForm(space, space, factory, assemblyOptions);
    DiscreteLinearOperatorPtr discreteDlp =
            dlp.assembleWeakForm(space, space, factory, assemblyOptions);
    DiscreteLinearOperatorPtr discreteId =
            id.assembleWeakForm(space, space, factory, assemblyOptions);

    arma::Mat<double> slpMatrix = discreteSlp->asMatrix();
    arma::Mat<double> dlpMatrix = discreteDlp->asMatrix();
    arma::Mat<double> idMatrix = discreteId->asMatrix();

    arma::Mat<double> lhs = slpMatrix;
    arma::Col<double> V(idMatrix.n_rows);
    V.fill(1.);
    arma::Col<double> rhs = (-0.5 * idMatrix + dlpMatrix) * V;
    solutionCoefficients = arma::solve(lhs, rhs);

    GridFunction<double> trace0(space, V);
    GridFunction<double> trace1(space, solutionCoefficients);

    trace1.exportToVtk(VtkWriter::VERTEX_DATA, "Neumann_data", "physical_test_neumann_data_vertex");
    // doesn't make very much sense -- just for testing
    trace1.exportToVtk(VtkWriter::CELL_DATA, "Neumann_data", "physical_test_neumann_data_cell");

    GridParameters params;
    params.topology = GridParameters::TETRAHEDRAL;

    std::cout << "Importing 3D grid" << std::endl;
    std::auto_ptr<Grid> volumeGrid = GridFactory::importGmshGrid(
                params,
                "spherical_shell_2194_nodes_volume.msh",
                true, // verbose
                false); // insertBoundarySegments

    std::auto_ptr<InterpolatedFunction<double> > slpOfTrace1 =
            slp.applyOffSurface(trace1, *volumeGrid, factory, evaluationOptions);
    std::auto_ptr<InterpolatedFunction<double> > dlpOfTrace0 =
            dlp.applyOffSurface(trace0, *volumeGrid, factory, evaluationOptions);
    slpOfTrace1->exportToVtk("slp_of_trace1", "physical_test_slp_of_trace1");
    dlpOfTrace0->exportToVtk("dlp_of_trace0", "physical_test_dlp_of_trace0");

    solutionCoefficients.print("Solution:\n");
    arma::Col<double> deviation = solutionCoefficients - (-1.);
    // % in Armadillo -> elementwise multiplication
    double stdDev = sqrt(arma::accu(deviation % deviation) /
                         solutionCoefficients.n_rows);
    std::cout << "Standard deviation: " << stdDev << std::endl;
}
