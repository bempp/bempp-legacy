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

#include "config_alugrid.hpp"
#include "config_trilinos.hpp"

#include "meshes.hpp"

#include "assembly/assembly_options.hpp"
#include "assembly/discrete_linear_operator.hpp"
#include "assembly/evaluation_options.hpp"
#include "assembly/grid_function.hpp"
#include "assembly/interpolated_function.hpp"
#include "assembly/linear_operator_superposition.hpp"
#include "assembly/ordinary_function.hpp"
#include "assembly/standard_local_assembler_factory_for_operators_on_surfaces.hpp"

#include "assembly/identity_operator.hpp"
#include "assembly/single_layer_potential_3d.hpp"
#include "assembly/double_layer_potential_3d.hpp"
#include "assembly/adjoint_double_layer_potential_3d.hpp"
#include "assembly/hypersingular_operator_3d.hpp"

#include "grid/geometry.hpp"
#include "grid/geometry_factory.hpp"
#include "grid/grid_factory.hpp"

#include "linalg/aca_preconditioner_factory.hpp"
#include "linalg/default_iterative_solver.hpp"
#include "linalg/default_direct_solver.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"

#include <armadillo>
#include <cmath>
#include <iostream>
#include <memory>

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
        //result(0) = sin(0.5 * point(0) * cos(0.25 * point(2))) * cos(point(1));
        result(0) = 1.;
    }
};


int main()
{
    // Load a predefined test grid
    const MeshVariant meshVariant = SPHERE_614;
    std::auto_ptr<Grid> grid = loadMesh(meshVariant);

    // Initialize the spaces

    PiecewiseLinearContinuousScalarSpace<double> HplusHalfSpace(*grid);
    PiecewiseConstantScalarSpace<double> HminusHalfSpace(*grid);

    HplusHalfSpace.assignDofs();
    HminusHalfSpace.assignDofs();

    // Define some default options.

    AssemblyOptions assemblyOptions;
    EvaluationOptions evaluationOptions;

    // We want to use ACA

    AcaOptions acaOptions; // Default parameters for ACA
    //assemblyOptions.switchToAca(acaOptions);

    // Define the standard integration factory

    StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double> factory;

    // We need the single layer, double layer, and the identity operator

    SingleLayerPotential3D<double> slp(HminusHalfSpace,HminusHalfSpace);
    DoubleLayerPotential3D<double> dlp(HminusHalfSpace,HplusHalfSpace);
    IdentityOperator<double> id(HminusHalfSpace,HplusHalfSpace);

    // Form the right-hand side sum

    LinearOperatorSuperposition<double> rhsOp = -0.5 * id + dlp;

    // Assemble the Operators

    slp.assemble(factory,assemblyOptions);
    rhsOp.assemble(factory,assemblyOptions);

    // We also want a grid function

    MyFunctor functor;
    OrdinaryFunction<MyFunctor> function(functor);

    GridFunction<double> u(HplusHalfSpace, function, factory, assemblyOptions);

    // Assemble the rhs

    std::cout << "Assemble rhs" << std::endl;

    GridFunction<double> rhs = rhsOp * u;

    // Initialize the iterative solver

    std::cout << "Initialize solver" << std::endl;

#ifdef WITH_TRILINOS
    DefaultIterativeSolver<double> solver(slp, rhs);
    solver.initializeSolver(defaultGmresParameterList(1e-5));
    solver.solve();
    std::cout << solver.getSolverMessage() << std::endl;
#else
    DefaultDirectSolver<double> solver(slp, rhs);
    solver.solve();
#endif
    // Extract the solution

    GridFunction<double> solFun = solver.getResult();

    // Write out as VTK

    solFun.exportToVtk(VtkWriter::CELL_DATA, "Neumann_data",
                       "physical_test_neumann_data");

    // Compare with exact solution

    arma::Col<double> solutionCoefficients = solFun.coefficients().asArmadilloVector();

    arma::Col<double> deviation = solutionCoefficients - (-1.);
    // % in Armadillo -> elementwise multiplication
    double stdDev = sqrt(arma::accu(deviation % deviation) /
                         solutionCoefficients.n_rows);
    std::cout << "Standard deviation: " << stdDev << std::endl;
}
