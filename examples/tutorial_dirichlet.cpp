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

#include "grid/grid_factory.hpp"

#include "assembly/assembly_options.hpp"
#include "assembly/discrete_linear_operator.hpp"

#include "assembly/identity_operator.hpp"
#include "assembly/single_layer_potential_3d.hpp"
#include "assembly/double_layer_potential_3d.hpp"

#include "fiber/standard_local_assembler_factory_for_operators_on_surfaces.hpp"

#include "grid/geometry.hpp"
#include "grid/geometry_factory.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"

using namespace Bempp;


/**
  Solve the Laplace equation

    \nabla^2 V = 0

  outside a sphere with radius 1, with the Dirichlet boundary condition

    V(x) = 1

  imposed on the surface of the sphere and the condition

    lim_{|x|\to\infty} V(x) = 0

  imposed at infinity. The analytical solution is

    V(x) = 1/|x|,

  where |x| is the distance from the centre of the sphere. Thus, the derivative
  of V on the surface in the (outer) normal direction is -1.  More details about this example are given
  in the tutorial.
  */
void tutorial_dirichlet1(){
	// Load the mesh
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    std::auto_ptr<Grid> grid = GridFactory::importGmshGrid(params, "sphere-614.msh", false, false);

    // Define the approximation spaces
    PiecewiseLinearContinuousScalarSpace<double> linearCtsSpace(*grid);
    PiecewiseConstantScalarSpace<double> constSpace(*grid);

    linearCtsSpace.assignDofs();
    constSpace.assignDofs();

    // Representations of the continuous operators
    SingleLayerPotential3D<double> slp;
    DoubleLayerPotential3D<double> dlp;
    IdentityOperator<double> id;

    // Initialise Fiber
    Fiber::AccuracyOptions accuracyOptions; // just use the default options
    Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory> factory(accuracyOptions);

    AssemblyOptions assemblyOptions; // again, use the default
    typedef std::auto_ptr<DiscreteLinearOperator<double> > DiscreteLinearOperatorPtr;
    // build the weak forms of the operators
    DiscreteLinearOperatorPtr discreteSlp =
			slp.assembleWeakForm(constSpace, linearCtsSpace, factory, assemblyOptions);
	DiscreteLinearOperatorPtr discreteDlp =
			dlp.assembleWeakForm(constSpace, linearCtsSpace, factory, assemblyOptions);
	DiscreteLinearOperatorPtr discreteId =
			id.assembleWeakForm(constSpace, linearCtsSpace, factory, assemblyOptions);

	// Assemble the linear system ...
	arma::Mat<double> M = -0.5 * discreteId->asMatrix() + discreteDlp->asMatrix();
	arma::Mat<double> V = discreteSlp->asMatrix();
    arma::Col<double> g = arma::ones(linearCtsSpace.globalDofCount(), 1);

	// and solve it
    arma::Col<double> b = arma::solve(V, M * g);
	std::cout<<"Neumann coefficients"<<b.t();

}

int main() {
	tutorial_dirichlet1();
}
