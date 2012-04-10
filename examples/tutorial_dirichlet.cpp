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
#include "assembly/discrete_scalar_valued_linear_operator.hpp"

#include "assembly/identity_operator.hpp"
#include "assembly/single_layer_potential_3d.hpp"
#include "assembly/double_layer_potential_3d.hpp"
#include "assembly/adjoint_double_layer_potential_3d.hpp"
#include "assembly/hypersingular_operator_3d.hpp"
#include "assembly/source_term.hpp"
#include "assembly/discrete_scalar_valued_source_term.hpp"

#include "fiber/standard_local_assembler_factory_for_operators_on_surfaces.hpp"
#include "fiber/ordinary_function.hpp"

#include "grid/geometry.hpp"
#include "grid/geometry_factory.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"

using namespace Bempp;

class MyFunctor
{
public:
    // Type of the function's values (e.g. float or std::complex<double>)
    typedef double ValueType;

    // Number of components of the function's argument
    static const int argumentDimension = 3;

    // Number of components of the function's result
    static const int resultDimension = 1;

    // Evaluate the function at the point "point" and store result in
    // the array "result"
    inline void evaluate(const arma::Col<ValueType>& point,
                  arma::Col<ValueType>& result) const {
//        result(0) = sin(0.5 * point(0) * cos(0.25 * point(2))) * cos(point(1));
        result(0) = 1.0;
    }
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
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;

//    std::auto_ptr<Grid> grid = GridFactory::importGmshGrid(params, "sphere-614.msh", false, false);
    std::auto_ptr<Grid> grid = GridFactory::importGmshGrid(params, "cube-12-reoriented.msh", false, false);

    PiecewiseLinearContinuousScalarSpace<double> space(*grid);
    space.assignDofs();

    AssemblyOptions assemblyOptions;
//    assemblyOptions.switchToDense();
//    assemblyOptions.switchToTbb();
//    assemblyOptions.setSingularIntegralCaching(AssemblyOptions::YES);

//    Fiber::OpenClOptions openClOptions;
    //    assemblyOptions.switchToOpenCl(openClOptions);

    Fiber::AccuracyOptions accuracyOptions; // default
    Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory>
            factory(accuracyOptions);

    typedef std::auto_ptr<LinearOperator<double> > LinearOperatorPtr;

    arma::Col<double> solution;
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

//	MyFunctor myfn;
//    Fiber::OrdinaryFunction<MyFunctor> function(myfn);
//    SourceTerm<double> sourceTerm;
//    std::auto_ptr<DiscreteScalarValuedSourceTerm<double> > result =
//            sourceTerm.assembleWeakForm(function, space, factory, assemblyOptions);

    arma::Col<double> g = arma::ones(idMatrix.n_cols);
    std::cout<<g;
    std::cout<<"ones"<<idMatrix * arma::ones(idMatrix.n_cols, 1);

    arma::Col<double> rhs = (-0.5 * idMatrix + dlpMatrix) * g;
	solution = arma::solve(lhs, rhs);
	std::cout<<"solution"<<solution;

//	arma::Col<double> deviation = solution - (-1.);
//    // % in Armadillo -> elementwise multiplication
//    double stdDev = sqrt(arma::accu(deviation % deviation) / solution.n_rows);
//    std::cout << "Standard deviation: " << stdDev << std::endl;
}
