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

#include "assembly/identity_operator.hpp"
#include "assembly/single_layer_potential_3d.hpp"
#include "assembly/double_layer_potential_3d.hpp"
#include "assembly/adjoint_double_layer_potential_3d.hpp"
#include "assembly/hypersingular_operator_3d.hpp"

#include "fiber/standard_local_assembler_factory_for_operators_on_surfaces.hpp"

#include "grid/geometry.hpp"
#include "grid/geometry_factory.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"

using namespace Bempp;

enum Formulation
{
    CBIE,
    HBIE // currently doesn't work correctly -- to investigate!
};

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
    const MeshVariant meshVariant = SPHERE_152;
    const Formulation formulation = CBIE;

    std::auto_ptr<Grid> grid = loadMesh(meshVariant);
    dumpElementList(grid.get());

    PiecewiseLinearContinuousScalarSpace<double> space(*grid);
    space.assignDofs();

    AssemblyOptions assemblyOptions;
    assemblyOptions.switchToDense();

    assemblyOptions.switchToTbb();
    assemblyOptions.setSingularIntegralCaching(AssemblyOptions::YES);

    Fiber::OpenClOptions openClOptions;
    //    assemblyOptions.switchToOpenCl(openClOptions);

    Fiber::AccuracyOptions accuracyOptions; // default
    Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory>
            factory(accuracyOptions);

    typedef std::auto_ptr<LinearOperator<double> > LinearOperatorPtr;

    arma::Col<double> solution;
    if (formulation == CBIE)
    {
        LinearOperatorPtr slp(new SingleLayerPotential3D<double>);
        LinearOperatorPtr dlp(new DoubleLayerPotential3D<double>);
        LinearOperatorPtr id(new IdentityOperator<double>);

        typedef DiscreteScalarValuedLinearOperator<double> DiscreteLinearOperator;
        typedef std::auto_ptr<DiscreteLinearOperator> DiscreteLinearOperatorPtr;
        DiscreteLinearOperatorPtr discreteSlp =
                slp->assembleWeakForm(space, space, factory, assemblyOptions);
        DiscreteLinearOperatorPtr discreteDlp =
                dlp->assembleWeakForm(space, space, factory, assemblyOptions);
        DiscreteLinearOperatorPtr discreteId =
                id->assembleWeakForm(space, space, factory, assemblyOptions);

        arma::Mat<double> slpMatrix = discreteSlp->asMatrix();
        arma::Mat<double> dlpMatrix = discreteDlp->asMatrix();
        arma::Mat<double> idMatrix = discreteId->asMatrix();

        arma::Mat<double> lhs = slpMatrix;
        arma::Col<double> V(idMatrix.n_rows);
        V.fill(1.);
        arma::Col<double> rhs = (-0.5 * idMatrix + dlpMatrix) * V;
        solution = arma::solve(lhs, rhs);
    }
    else if (formulation == HBIE)
    {
        LinearOperatorPtr adlp(new AdjointDoubleLayerPotential3D<double>);
        LinearOperatorPtr hso(new HypersingularOperator3D<double>);
        LinearOperatorPtr id(new IdentityOperator<double>);

        typedef DiscreteScalarValuedLinearOperator<double> DiscreteLinearOperator;
        typedef std::auto_ptr<DiscreteLinearOperator> DiscreteLinearOperatorPtr;
        DiscreteLinearOperatorPtr discreteAdlp =
                adlp->assembleWeakForm(space, space, factory, assemblyOptions);
        DiscreteLinearOperatorPtr discreteHso =
                hso->assembleWeakForm(space, space, factory, assemblyOptions);
        DiscreteLinearOperatorPtr discreteId =
                id->assembleWeakForm(space, space, factory, assemblyOptions);

        arma::Mat<double> adlpMatrix = discreteAdlp->asMatrix();
        arma::Mat<double> hsoMatrix = discreteHso->asMatrix();
        arma::Mat<double> idMatrix = discreteId->asMatrix();

        arma::Mat<double> lhs = 0.5 * idMatrix - adlpMatrix;
        arma::Col<double> V(idMatrix.n_rows);
        V.fill(1.);
        arma::Col<double> rhs = -hsoMatrix * V;

        hsoMatrix.print("Hypersingular operator");
        rhs.print("Right-hand side");
        solution = arma::solve(lhs, rhs);
    }
    else
        throw std::runtime_error("Invalid formulation");

    solution.print("Solution:\n");
    arma::Col<double> deviation = solution - (-1.);
    // % in Armadillo -> elementwise multiplication
    double stdDev = sqrt(arma::accu(deviation % deviation) / solution.n_rows);
    std::cout << "Standard deviation: " << stdDev << std::endl;
}
