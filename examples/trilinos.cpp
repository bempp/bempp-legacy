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

#include "config_trilinos.hpp"
#include "config_ahmed.hpp"

#ifdef WITH_TRILINOS
#ifdef WITH_AHMED

#include "meshes.hpp"

#include "assembly/aca_approximate_lu_inverse.hpp"
#include "assembly/assembly_options.hpp"
#include "assembly/discrete_aca_linear_operator.hpp"
#include "assembly/discrete_linear_operator.hpp"
#include "assembly/vector.hpp"
#include "space/piecewise_constant_scalar_space.hpp"

#include "assembly/identity_operator.hpp"
#include "assembly/double_layer_potential_3d.hpp"
#include "assembly/single_layer_potential_3d.hpp"
#include "assembly/linear_operator_superposition.hpp"
#include "linalg/aca_preconditioner_factory.hpp"
#include "linalg/default_iterative_solver.hpp"

#include "common/auto_timer.hpp"

#include "fiber/ordinary_function.hpp"
#include "fiber/standard_local_assembler_factory_for_operators_on_surfaces.hpp"

#include "grid/geometry.hpp"
#include "grid/geometry_factory.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"


#include <armadillo>
#include <cmath>
#include <iostream>
#include <memory>
#include <typeinfo>

#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Thyra_AztecOOLinearOpWithSolveFactory.hpp>
#include <Thyra_BelosLinearOpWithSolveFactory_decl.hpp>
#include <Thyra_DetachedVectorView.hpp>
#include <Thyra_LinearOpTester.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>

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
        result(0) = sin(0.5 * point(0) * cos(0.25 * point(2))) * cos(point(1));
    }
};

/**
    This script demonstrates how to use a solver from the Trilinos.Belos library
    to solve a system of equations generated with BEM++.
  */
int main()
{
    // Reference-counted pointer from Trilinos
    using Teuchos::RCP;

    // OPTIONS CONTROLLING THE EXECUTION OF THE SCRIPT
    // Mesh to use
    const MeshVariant meshVariant = SPHERE_152; // SPHERE_2590;
    // If true, test whether the discrete linear operator is implemented correctly
    const bool testLinearProperties = false;
    // If true, compare results obtained with the direct and iterative solvers
    const bool compareArmadilloAndBelos = true;
    // If true, assembly in ACA mode, otherwise in dense mode
    const bool useAca = true;

    std::auto_ptr<Grid> grid = loadMesh(meshVariant);

    //PiecewiseLinearContinuousScalarSpace<double> space(*grid);
    PiecewiseConstantScalarSpace<double> space(*grid);

    space.assignDofs();

    AssemblyOptions assemblyOptions;

    AcaOptions acaOptions;
    acaOptions.eps=1E-3;
    acaOptions.eta=1.2;
/*
    acaOptions.eps = 1e-4;
    acaOptions.maximumRank = 10000;
    acaOptions.minimumBlockSize = 32;
    acaOptions.eta = 0.8;
    acaOptions.recompress = true;
*/

    assemblyOptions.switchToAca(acaOptions);
    assemblyOptions.switchToTbb();
    assemblyOptions.setSingularIntegralCaching(AssemblyOptions::YES);

    Fiber::AccuracyOptions accuracyOptions; // default
    //accuracyOptions.regular.mode=accuracyOptions.regular.EXACT_ORDER;
    //accuracyOptions.regular.order=8;

    Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory>
            factory(accuracyOptions);

    // Create formal linear operators

    typedef std::auto_ptr<ElementaryLinearOperator<double> > LinearOperatorPtr;   
    LinearOperatorPtr halfId(new IdentityOperator<double>);
    halfId->scale(0.5);
    //LinearOperatorPtr potential(new DoubleLayerPotential3D<double>);
    SingleLayerPotential3D<double> lhsOperator;

    //LinearOperatorSuperposition<double> lhsOperator(halfId, potential);

    // Generate discrete linear operator

    typedef DiscreteLinearOperator<double> DiscreteLinearOperator;
    typedef std::auto_ptr<DiscreteLinearOperator> DiscreteLinearOperatorPtr;
    DiscreteLinearOperatorPtr discreteLhs =
            lhsOperator.assembleWeakForm(space, space, factory, assemblyOptions);

    // Create a right-hand side

    MyFunctor functor;
    Fiber::OrdinaryFunction<MyFunctor> function(functor);

    typedef std::auto_ptr<Vector<double> > VectorPtr;

    arma::Col<double> rhsVector(space.globalDofCount());
    rhsVector.fill(1.);
    VectorPtr discreteRhs(new Vector<double>(rhsVector));

    DefaultIterativeSolver<double> iterativeSolver(*discreteLhs,*discreteRhs);

    // It is also possible to initialize with a vector of right-hand sides for block-solves, e.g.
    // std::vector<Vector<double>* > srcTerms;
    // srcTerms.push_back(&(*discreteRhs));
    // DefaultIterativeSolver<double> iterativeSolver(*discreteLhs,srcTerms);


    //iterativeSolver.addPreconditioner(AcaPreconditionerFactory<double>::AcaOperatorToPreconditioner(*discreteLhs,0.1));
    iterativeSolver.initializeSolver(defaultGmresParameterList(1E-8));
    //iterativeSolver.initializeSolver(Teuchos::getParametersFromXmlFile("trilinos-belos.xml"));
    iterativeSolver.solve();
    std::cout << iterativeSolver.getSolverMessage() << std::endl;

    // Compare the solution produced by Belos with one obtained with Armadillo
    if (compareArmadilloAndBelos)
    {
        // Solve the system using Armadillo
        arma::Mat<double> lhsMatrix = discreteLhs->asMatrix();
        arma::Col<double> armaSolution = arma::solve(lhsMatrix, rhsVector);

        arma::Mat<double> belosSolution=iterativeSolver.getResult();

        std::cout << belosSolution.n_rows << std::endl;
        std::cout << belosSolution.n_cols << std::endl;


        double diff = arma::norm(armaSolution - belosSolution, 2);
        std::cout << "Norm of difference between Armadillo and Belos solutions: "
                  << diff << std::endl;
    }
}


#else

#include <iostream>
int main(void){
    std::cout << "This example requires Trilinos and AHMED " << std::endl;
}
#endif // WITH_AHMED
#else
#include <iostream>
int main(void){
    std::cout << "This example requires Trilinos and AHMED " << std::endl;
}
#endif // WITH_TRILINOS
