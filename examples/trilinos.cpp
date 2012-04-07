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

#include "../lib/common/config_trilinos.hpp"
#include "../lib/common/config_ahmed.hpp"

#ifdef WITH_TRILINOS

#include "meshes.hpp"

#include "assembly/aca_approximate_lu_inverse.hpp"
#include "assembly/assembly_options.hpp"
#include "assembly/discrete_aca_scalar_valued_linear_operator.hpp"
#include "assembly/discrete_scalar_valued_linear_operator.hpp"
#include "assembly/discrete_scalar_valued_source_term.hpp"
#include "assembly/source_term.hpp"

#include "assembly/identity_operator.hpp"
#include "assembly/double_layer_potential_3d.hpp"
#include "assembly/linear_operator_superposition.hpp"

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
    const MeshVariant meshVariant = CUBE_384;
    // If true, test whether the discrete linear operator is implemented correctly
    const bool testLinearProperties = false;
    // If true, compare results obtained with the direct and iterative solvers
    const bool compareArmadilloAndBelos = true;
    // If true, assembly in ACA mode, otherwise in dense mode
    const bool useAca = false;

    std::auto_ptr<Grid> grid = loadMesh(meshVariant);
    dumpElementList(grid.get());

    PiecewiseLinearContinuousScalarSpace<double> space(*grid);
    space.assignDofs();

    AssemblyOptions assemblyOptions;

    AcaOptions acaOptions;
    acaOptions.eps = 1e-4;
    acaOptions.maximumRank = 10000;
    acaOptions.minimumBlockSize = 2;
    acaOptions.eta = 0.8;
    acaOptions.recompress = true;
    if (useAca)
        assemblyOptions.switchToAca(acaOptions);
    else
        assemblyOptions.switchToDense();

    assemblyOptions.switchToTbb();
    assemblyOptions.setSingularIntegralCaching(AssemblyOptions::YES);

    Fiber::AccuracyOptions accuracyOptions; // default
    Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory>
            factory(accuracyOptions);

    // Create formal linear operators

    typedef std::auto_ptr<ElementaryLinearOperator<double> > LinearOperatorPtr;   
    LinearOperatorPtr halfId(new IdentityOperator<double>);
    halfId->scale(0.5);
    LinearOperatorPtr potential(new DoubleLayerPotential3D<double>);

    LinearOperatorSuperposition<double> lhsOperator(halfId, potential);

    // Generate discrete linear operator

    typedef DiscreteScalarValuedLinearOperator<double> DiscreteLinearOperator;
    typedef std::auto_ptr<DiscreteLinearOperator> DiscreteLinearOperatorPtr;
    DiscreteLinearOperatorPtr discreteLhs =
            lhsOperator.assembleWeakForm(space, space, factory, assemblyOptions);

    // Create a right-hand side

    MyFunctor functor;
    Fiber::OrdinaryFunction<MyFunctor> function(functor);

    SourceTerm<double> sourceTerm;
    typedef DiscreteScalarValuedSourceTerm<double> DiscreteSourceTerm;
    typedef std::auto_ptr<DiscreteSourceTerm> DiscreteSourceTermPtr;

    DiscreteSourceTermPtr discreteRhs =
            sourceTerm.assembleWeakForm(function, space, factory, assemblyOptions);

    // Test the discrete linear operator

    // Output stream
    RCP<Teuchos::FancyOStream> out =
            Teuchos::VerboseObjectBase::getDefaultOStream();

    if (testLinearProperties)
    {
        const double tolerance = 1e-8;
        Thyra::LinearOpTester<double> linearOpTester;
        linearOpTester.enable_all_tests(false);
        linearOpTester.check_linear_properties(true);
        linearOpTester.set_all_error_tol(tolerance);
        linearOpTester.set_all_warning_tol(1e-2 * tolerance);
        linearOpTester.show_all_tests(true);
        linearOpTester.check(*discreteLhs, out.ptr());
    }

    
    std::cout << "Operator tested " << std::endl;
    
    //std::cout << discreteLhs->asMatrix() << std::endl;
    
    // Solve the resulting system

    // Create a factory for invertible operators
    Thyra::BelosLinearOpWithSolveFactory<double> invertibleOpFactory;
    RCP<Teuchos::ParameterList> paramList =
            Teuchos::getParametersFromXmlFile("trilinos-belos.xml");
    invertibleOpFactory.setParameterList(paramList);
    invertibleOpFactory.setOStream(out);
    invertibleOpFactory.setVerbLevel(Teuchos::VERB_HIGH);

    // Wrap the discrete operator in a Trilinos reference-counted pointer
    // (non-owning -- destruction is managed by the original auto_ptr)
    RCP<const Thyra::LinearOpBase<double> > trilinosDiscreteLhs(
                discreteLhs.get(), false /*don't own*/);

    // Upgrade the discrete linear operator to an object supporting the solve operation

    RCP<Thyra::LinearOpWithSolveBase<double> > invertibleDiscreteLhs;
    
    std::cout << "Invertible Operator created" << std::endl;

#ifdef WITH_AHMED
    if (assemblyOptions.operatorRepresentation() == AssemblyOptions::ACA)
    {
        try {
        // (a) variant with preconditioner
        AutoTimer precTimer("Preconditioner building took ");
        const double delta = 0.1; // LU approximation accuracy
        // We know that we've built an ACA representation of the operator,
        // so this dynamic cast will succeed
        std::cout << "Dynamic Cast" << std::endl;
        const DiscreteAcaScalarValuedLinearOperator<double>& discreteAcaLhs = 
        DiscreteAcaScalarValuedLinearOperator<double>::castToAca(*discreteLhs);
        // Construct the linear operator to act as a preconditioner
        std::cout << "Start creating preconditioner" << std::endl;
        RCP<const Thyra::LinearOpBase<double> > precOp(
                    new AcaApproximateLuInverse<double>(discreteAcaLhs, delta));
        // and wrap it in a PreconditionerBase object
        // (the static cast is there because unspecifiedPrec() returns
        // a ref-counted pointer to a subclass of PreconditionerBase
        std::cout << "Created approximate inverse" << std::endl;
        RCP<const Thyra::PreconditionerBase<double> > preconditioner =
                Teuchos::rcp_static_cast<const Thyra::PreconditionerBase<double> >(
                    Thyra::unspecifiedPrec(precOp));
        // Now create a discrete linear operator with a solve() member function
        invertibleDiscreteLhs = invertibleOpFactory.createOp();
        Thyra::initializePreconditionedOp(
                    invertibleOpFactory,
                    trilinosDiscreteLhs,
                    preconditioner,
                    invertibleDiscreteLhs.ptr());
        }
        catch( const std::bad_cast &e ) {
            std::cout << e.what() << std::endl;
        }
    }
    else
    {
#endif // WITH_AHMED
        // (b) variant without preconditioner
        // Create a discrete linear operator with a solve() member function
        invertibleDiscreteLhs =
                Thyra::linearOpWithSolve(invertibleOpFactory, trilinosDiscreteLhs);
#ifdef WITH_AHMED
    }
#endif // WITH_AHMED

    // Create the solution vector
    RCP<Thyra::VectorBase<double> > solution =
            Thyra::createMember(invertibleDiscreteLhs->domain());
    // and fill it with zeros
    Thyra::assign(solution.ptr(), 0.);
    // Solve the system using default solve criteria
    Thyra::SolveStatus<double> status =
            Thyra::solve<double>(*invertibleDiscreteLhs, Thyra::NOTRANS,
                                 *discreteRhs, solution.ptr());
    std::cout << "\nSolve status:\n" << status;
    std::cout << "Solution: " << *solution << "\n";

    // Compare the solution produced by Belos with one obtained with Armadillo
    if (compareArmadilloAndBelos)
    {
        // Solve the system using Armadillo
        arma::Mat<double> lhsMatrix = discreteLhs->asMatrix();
        arma::Col<double> rhsVector = discreteRhs->asVector();
        arma::Col<double> armaSolution = arma::solve(lhsMatrix, rhsVector);

        // Wrap Belos's solution in an Armadillo vector
        Thyra::ConstDetachedVectorView<double> belosView(solution);
        arma::Col<double> belosSolution(belosView.values(), belosView.subDim());

        double diff = arma::norm(armaSolution - belosSolution, 2);
        std::cout << "Norm of difference between Armadillo and Belos solutions: "
                  << diff << std::endl;
    }
}


#else

#include <iostream>
int main(void){
    std::cout << "This example requires trilinos" << std::endl;
}
#endif // WITH_TRILINOS
