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
#include "assembly/standard_local_assembler_factory_for_operators_on_surfaces.hpp"

#include "assembly/identity_operator.hpp"
#include "assembly/laplace_3d_single_layer_potential.hpp"
#include "assembly/laplace_3d_double_layer_potential.hpp"

#include "common/scalar_traits.hpp"

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

typedef double BFT; // basis function type
typedef double RT; // result type (type used to represent discrete operators)

class MyFunctor
{
public:
    // Type of the function's values (e.g. float or std::complex<double>)
    typedef RT ValueType;
    // Type of coordinates (must be the "real part" of ValueType)
    typedef ScalarTraits<RT>::RealType CoordinateType;

    // Number of components of the function's argument
    static const int argumentDimension = 3;

    // Number of components of the function's result
    static const int resultDimension = 1;

    // Evaluate the function at the point "point" and store result in
    // the array "result"
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        result(0) = 1.;
    }
};


int main(int argc, char* argv[])
{
    // Load mesh

    if (argc != 2) {
        std::cout << "Solve a Dirichlet problem for the Laplace equation.\n"
                     "Usage: " << argv[0] << " <mesh_file>" << std::endl;
        return 1;
    }
    std::auto_ptr<Grid> grid = loadTriangularMeshFromFile(argv[1]);

    // Initialize the spaces

    PiecewiseLinearContinuousScalarSpace<BFT> HplusHalfSpace(*grid);
    PiecewiseConstantScalarSpace<BFT> HminusHalfSpace(*grid);

    HplusHalfSpace.assignDofs();
    HminusHalfSpace.assignDofs();

    // Define some default options.

    AssemblyOptions assemblyOptions;

    // We want to use ACA

    AcaOptions acaOptions; // Default parameters for ACA
    acaOptions.eps = 1e-5;
    assemblyOptions.switchToAca(acaOptions);

    // Define the standard integration factory

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.orderIncrement = 1;
    StandardLocalAssemblerFactoryForOperatorsOnSurfaces<BFT, RT>
            factory(accuracyOptions);

    // We need the single layer, double layer, and the identity operator

    Laplace3dSingleLayerPotential<BFT, RT> slp(HminusHalfSpace, HminusHalfSpace);
    Laplace3dDoubleLayerPotential<BFT, RT> dlp(HminusHalfSpace, HplusHalfSpace);
    IdentityOperator<BFT, RT> id(HminusHalfSpace, HplusHalfSpace);

    // Form the right-hand side sum

    LinearOperatorSuperposition<BFT, RT> rhsOp = -0.5 * id + dlp;

    // Assemble the operators

    slp.assemble(factory, assemblyOptions);
    rhsOp.assemble(factory, assemblyOptions);

    // We also want a grid function

    GridFunction<BFT, RT> u = surfaceNormalIndependentGridFunction(
                HplusHalfSpace, MyFunctor(), factory, assemblyOptions);

    // Assemble the rhs

    std::cout << "Assemble rhs" << std::endl;

    GridFunction<BFT, RT> rhs = rhsOp * u;

    // Initialize the solver

    std::cout << "Initialize solver" << std::endl;

    // WARNING: The default iterative solver does not support
    // complex-valued operators, hence if you typedef RT to std::complex<...>,
    // you need to use the direct solver.
#ifdef WITH_TRILINOS
    DefaultIterativeSolver<BFT, RT> solver(slp, rhs);
    solver.initializeSolver(defaultGmresParameterList(1e-5));
    solver.solve();
    std::cout << solver.getSolverMessage() << std::endl;
#else
    DefaultDirectSolver<BFT, RT> solver(slp, rhs);
    solver.solve();
#endif

    // Extract the solution

    GridFunction<BFT, RT> solFun = solver.getResult();

    // Write out as VTK

    solFun.exportToVtk(VtkWriter::CELL_DATA, "Neumann_data",
                       "calculated_neumann_data_cell");
    solFun.exportToVtk(VtkWriter::VERTEX_DATA, "Neumann_data",
                       "calculated_neumann_data_vertex");

//    // Uncomment the block below if you are solving the problem on a sphere and
//    // you want to compare the numerical and analytical solution.

//    arma::Col<RT> solutionCoefficients =
//            solFun.coefficients().asArmadilloVector();
//    std::cout << solutionCoefficients << std::endl;

//    arma::Col<RT> deviation = solutionCoefficients - static_cast<RT>(-1.);
//    // % in Armadillo -> elementwise multiplication
//    RT stdDev = sqrt(arma::accu(deviation % deviation) /
//                     static_cast<RT>(solutionCoefficients.n_rows));
//    std::cout << "Standard deviation: " << stdDev << std::endl;
}
