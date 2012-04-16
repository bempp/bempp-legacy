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

#include "assembly/identity_operator.hpp"
#include "assembly/single_layer_potential_3d.hpp"
#include "assembly/double_layer_potential_3d.hpp"
#include "assembly/adjoint_double_layer_potential_3d.hpp"
#include "assembly/hypersingular_operator_3d.hpp"

#include "fiber/standard_local_assembler_factory_for_operators_on_surfaces.hpp"
#include "fiber/ordinary_function.hpp"

#include "grid/geometry.hpp"
#include "grid/geometry_factory.hpp"
#include "grid/grid_factory.hpp"

#include "linalg/aca_preconditioner_factory.hpp"
#include "linalg/default_iterative_solver.hpp"

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
    assemblyOptions.switchToAca(acaOptions);

    // Define the standard integration factory



    Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory>
            factory;

    // We need the single layer, double layer, and the identity operator

    SingleLayerPotential3D<double> slp(HminusHalfSpace,HminusHalfSpace);
    DoubleLayerPotential3D<double> dlp(HminusHalfSpace,HplusHalfSpace);
    IdentityOperator<double> id(HminusHalfSpace,HplusHalfSpace);

    // Form the right-hand side sum

    LinearOperatorSuperposition<double> rhsOp=-.5*id+dlp;

    // Assemble the Operators

    slp.assemble(factory,assemblyOptions);
    rhsOp.assemble(factory,assemblyOptions);

    // We also want a grid function

    MyFunctor functor;
    Fiber::OrdinaryFunction<MyFunctor> function(functor);

    GridFunction<double> u(HplusHalfSpace, function, factory, assemblyOptions);


    // Assemble the rhs

    std::cout << "Assemble rhs" << std::endl;

    GridFunction<double> rhs=rhsOp.apply(u);

    // Initialize the iterative solver

    std::cout << "Initialize iterative solver" << std::endl;

    DefaultIterativeSolver<double> iterativeSolver(slp, rhs);
    iterativeSolver.initializeSolver(defaultGmresParameterList(1e-5));
    iterativeSolver.solve();
    std::cout << iterativeSolver.getSolverMessage() << std::endl;

    // Extract the solution

    GridFunction<double> solFun=iterativeSolver.getResult();

    // Write out as VTK

    solFun.exportToVtk(VtkWriter::CELL_DATA, "Neumann_data", "physical_test_neumann_data_vertex");

    // Compare with exact solution

    arma::Col<double> solutionCoefficients=solFun.coefficients().asArmadilloVector();

    arma::Col<double> deviation = solutionCoefficients - (-1.);
    // % in Armadillo -> elementwise multiplication
    double stdDev = sqrt(arma::accu(deviation % deviation) /
                         solutionCoefficients.n_rows);
    std::cout << "Standard deviation: " << stdDev << std::endl;





    // Finally discretize the operators

//    DiscreteLinearOperatorPtr discreteSlp =
//            slp.assembleWeakForm(factory, assemblyOptions);
//    DiscreteLinearOperatorPtr discreteRhsOp =
//            rhsOp.assembleWeakForm(factory, assemblyOptions);

//    slp.assemble(factory,assemblyOptions);

    // Now construct the right hand side

//    // Initialize the analytic rhs function

//    MyFunctor functor;
//    Fiber::OrdinaryFunction<MyFunctor> function(functor);

//    // A GridFunction is a discrete representation

//    GridFunction<double> trace0(HplusHalfSpace, function, factory, assemblyOptions);



//    // Extract the coefficient vector in the given basis

//    arma::Col<double> V = trace0.coefficients().asArmadilloVector();

//    std::cout << "Constructing RHS" << std::endl;

//    // Form (I+K)*b

//    arma::Col<double> rhs(HminusHalfSpace.globalDofCount());
//    discreteRhsOp->apply(NO_TRANSPOSE,V,rhs,1.,0.);

//    //discreteId->apply(NO_TRANSPOSE, V, rhs, -0.5, 0.);
//    //discreteDlp->apply(NO_TRANSPOSE, V, rhs, 1., 1.);

//    // Store in a discrete rhs vector

//    typedef std::auto_ptr<Vector<double> > VectorPtr;
//    VectorPtr discreteRhsVec(new Vector<double>(rhs));
//    GridFunction<double> rhs_grid(HminusHalfSpace,*discreteRhsVec);


//    // Initialize the iterative solver

//    DefaultIterativeSolver<double> iterativeSolver(slp, rhs_grid);

//    // We want a preconditioned solve

////    std::cout << "Constructing preconditioner" << std::endl;
////
////    iterativeSolver.addPreconditioner(
////                AcaPreconditionerFactory<double>::
////                AcaOperatorToPreconditioner(*discreteSlp, 0.1));

//    // Initialize the solver with parameters for GMRES

//    iterativeSolver.initializeSolver(defaultGmresParameterList(1e-8));
//    std::cout << "Solving" << std::endl;

//    // Solve and extract solution

//    iterativeSolver.solve();
//    std::cout << iterativeSolver.getSolverMessage() << std::endl;

//    arma::Mat<double> solutionCoefficients = iterativeSolver.getResult().coefficients().asArmadilloVector();

//    // We store the result into a Grid function

//    GridFunction<double> trace1(HminusHalfSpace, solutionCoefficients.col(0));

//    // Export the results to VTK

//    std::cout << "Exporting results to VTK" << std::endl;

//    trace1.exportToVtk(VtkWriter::VERTEX_DATA, "Neumann_data", "physical_test_neumann_data_vertex");


    // Compare with exact solution

//    arma::Col<double> deviation = solutionCoefficients - (-1.);
//    // % in Armadillo -> elementwise multiplication
//    double stdDev = sqrt(arma::accu(deviation % deviation) /
//                         solutionCoefficients.n_rows);
//    std::cout << "Standard deviation: " << stdDev << std::endl;
}
