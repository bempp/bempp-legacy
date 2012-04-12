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
        // result(0) = sin(0.5 * point(0) * cos(0.25 * point(2))) * cos(point(1));
        result(0) = 1.;
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
    const MeshVariant meshVariant = SPHERE_644;

    std::auto_ptr<Grid> grid = loadMesh(meshVariant);
    // dumpElementList(grid.get());

    PiecewiseLinearContinuousScalarSpace<double> HplusHalfSpace(*grid);
    PiecewiseConstantScalarSpace<double> HminusHalfSpace(*grid);
    HplusHalfSpace.assignDofs();
    HminusHalfSpace.assignDofs();

    AssemblyOptions assemblyOptions;
    AcaOptions acaOptions; // default
    assemblyOptions.switchToTbb();
    assemblyOptions.switchToAca(acaOptions);
    EvaluationOptions evaluationOptions;

    Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory>
            factory;

    typedef std::auto_ptr<LinearOperator<double> > LinearOperatorPtr;

    SingleLayerPotential3D<double> slp;
    DoubleLayerPotential3D<double> dlp;
    IdentityOperator<double> id;

    typedef DiscreteLinearOperator<double> DiscreteLinearOperator;
    typedef std::auto_ptr<DiscreteLinearOperator> DiscreteLinearOperatorPtr;
    DiscreteLinearOperatorPtr discreteSlp =
            slp.assembleWeakForm(HminusHalfSpace, HminusHalfSpace, factory, assemblyOptions);
    DiscreteLinearOperatorPtr discreteDlp =
            dlp.assembleWeakForm(HminusHalfSpace, HplusHalfSpace, factory, assemblyOptions);
#ifdef WITH_TRILINOS
    assemblyOptions.switchToSparse();
#endif
    DiscreteLinearOperatorPtr discreteId =
            id.assembleWeakForm(HminusHalfSpace, HplusHalfSpace, factory, assemblyOptions);

    MyFunctor functor;
    Fiber::OrdinaryFunction<MyFunctor> function(functor);

    GridFunction<double> trace0(HplusHalfSpace, function, factory, assemblyOptions);
    arma::Col<double> V = trace0.coefficients().asArmadilloVector();

    std::cout << "Constructing RHS" << std::endl;

    arma::Col<double> rhs(HminusHalfSpace.globalDofCount());
    discreteId->apply(NO_TRANSPOSE, V, rhs, -0.5, 0.);
    discreteDlp->apply(NO_TRANSPOSE, V, rhs, 1., 1.);

    typedef std::auto_ptr<Vector<double> > VectorPtr;
    VectorPtr discreteRhs(new Vector<double>(rhs));

#ifdef WITH_TRILINOS
    DefaultIterativeSolver<double> iterativeSolver(*discreteSlp, *discreteRhs);

    std::cout << "Constructing preconditioner" << std::endl;

    iterativeSolver.addPreconditioner(
                AcaPreconditionerFactory<double>::
                AcaOperatorToPreconditioner(*discreteSlp, 0.1));
    iterativeSolver.initializeSolver(defaultGmresParameterList(1e-5));
    std::cout << "Solving" << std::endl;
    iterativeSolver.solve();
    std::cout << iterativeSolver.getSolverMessage() << std::endl;

    arma::Mat<double> solutionCoefficients = iterativeSolver.getResult();
#else
    arma::Mat<double> solutionCoefficients = arma::solve(discreteSlp->asMatrix(), rhs);
#endif
    GridFunction<double> trace1(HminusHalfSpace, solutionCoefficients.col(0));

    std::cout << "Exporting results to VTK" << std::endl;

    trace1.exportToVtk(VtkWriter::VERTEX_DATA, "Neumann_data", "physical_test_neumann_data_vertex");
    // doesn't make very much sense -- just for testing
    trace1.exportToVtk(VtkWriter::CELL_DATA, "Neumann_data", "physical_test_neumann_data_cell");

// Plot field in volume
//    GridParameters params;
//    params.topology = GridParameters::TETRAHEDRAL;

//#ifdef WITH_ALUGRID
//    std::cout << "Importing 3D grid" << std::endl;
//    std::auto_ptr<Grid> volumeGrid = GridFactory::importGmshGrid(
//                params,
//                "spherical_shell_2194_nodes_volume.msh",
//                true, // verbose
//                false); // insertBoundarySegments

//    std::auto_ptr<InterpolatedFunction<double> > slpOfTrace1 =
//            slp.applyOffSurface(trace1, *volumeGrid, factory, evaluationOptions);
//    std::auto_ptr<InterpolatedFunction<double> > dlpOfTrace0 =
//            dlp.applyOffSurface(trace0, *volumeGrid, factory, evaluationOptions);
//    slpOfTrace1->exportToVtk("slp_of_trace1", "physical_test_slp_of_trace1");
//    dlpOfTrace0->exportToVtk("dlp_of_trace0", "physical_test_dlp_of_trace0");
//#endif

//    solutionCoefficients.print("Solution:\n");
    arma::Col<double> deviation = solutionCoefficients - (-1.);
    // % in Armadillo -> elementwise multiplication
    double stdDev = sqrt(arma::accu(deviation % deviation) /
                         solutionCoefficients.n_rows);
    std::cout << "Standard deviation: " << stdDev << std::endl;
}
