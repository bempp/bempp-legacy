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

#include "bempp/common/config_trilinos.hpp"

#include "../examples/meshes.hpp"

#include "assembly/assembly_options.hpp"
#include "assembly/blocked_boundary_operator.hpp"
#include "assembly/blocked_operator_structure.hpp"
#include "assembly/context.hpp"
#include "assembly/discrete_boundary_operator.hpp"
#include "assembly/grid_function.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "assembly/surface_normal_independent_function.hpp"

#include "assembly/identity_operator.hpp"
#include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_double_layer_boundary_operator.hpp"

#include "common/scalar_traits.hpp"

#include "linalg/aca_preconditioner_factory.hpp"
#include "linalg/default_iterative_solver.hpp"
#include "linalg/default_direct_solver.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"

#include <sys/wait.h>

using namespace Bempp;

typedef double BFT; // basis function type
typedef std::complex<double> RT; // result type (type used to represent discrete operators)

class NullFunctor
{
public:
    // Type of the function's values (e.g. float or std::complex<double>)
    typedef RT ValueType;
    // Type of coordinates (must be the "real part" of ValueType)
    typedef ScalarTraits<RT>::RealType CoordinateType;

    // Number of components of the function's argument
    int argumentDimension() const { return 3; }
    // Number of components of the function's result
    int resultDimension() const { return 1; }

    // Evaluate the function at the point "point" and store result in
    // the array "result"
    // Note that the source is located inside the inner sphere
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        result.fill(0.);
    }
};

class MyFunctor
{
public:
    // Type of the function's values (e.g. float or std::complex<double>)
    typedef RT ValueType;
    // Type of coordinates (must be the "real part" of ValueType)
    typedef ScalarTraits<RT>::RealType CoordinateType;

    // Number of components of the function's argument
    int argumentDimension() const { return 3; }
    // Number of components of the function's result
    int resultDimension() const { return 1; }

    MyFunctor(RT waveNumber) : m_waveNumber(waveNumber) {}

    // Evaluate the function at the point "point" and store result in
    // the array "result"
    // Note that the source is located inside the inner sphere
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        std::vector<CoordinateType> src(3);
        src[0] = 0; src[1] = -sqrt(60); src[2] = sqrt(60);
        CoordinateType dst = 0.0;
        for (int i = 0; i < 3; i++) {
            CoordinateType d = src[i] - point[i];
            dst += d*d;
        }
        dst = sqrt(dst);
        result(0) = exp (-m_waveNumber*dst);
	result(0) = 1;
    }

private:
    RT m_waveNumber;
};

template <typename T>
inline shared_ptr<T> SHARED(T& t)
{
    return make_shared_from_ref(t);
}

double A_Keijzer (double n)
{
    double th = asin(1.0/n);
    double costh = fabs(cos(th));
    double R0 = ((n-1.0)*(n-1.0)) / ((n+1.0)*(n+1.0));
    return (2.0/(1.0-R0) - 1.0 + costh*costh*costh) / (1.0 - costh*costh);
}

std::auto_ptr<Grid> CreateSphere (double rad, double elsize)
{
    // Create a sphere mesh of radius rad and mean element size elsize
    // Calls external gmsh executable which must be on the path

    const char *gmesh_def_name = "meshes/sphere.txt";
    const char *gmesh_geo_name = "meshes/sphere.geo";
    const char *gmesh_msh_name = "meshes/sphere.msh";
    char cbuf[512];

    std::ifstream ifs(gmesh_def_name);
    std::ofstream ofs(gmesh_geo_name);
    ofs << "rad = " << rad << ";\n";
    ofs << "lc = " << elsize << ";\n";
    while (ifs.getline (cbuf, 512))
        ofs << cbuf << std::endl;
    ofs.close();

    pid_t pID = vfork();
    if (pID == 0) {
        execlp ("gmsh", "gmsh", "-2", "meshes/sphere.geo", NULL);
    }
    int gmshExitStatus;
    waitpid (pID, &gmshExitStatus, 0);
    std::auto_ptr<Grid> grid = loadTriangularMeshFromFile(gmesh_msh_name);
    remove (gmesh_geo_name);
    remove (gmesh_msh_name);
    return grid;
}

int main(int argc, char* argv[])
{
    // Physical parameters, general
    const BFT c0 = 0.3;      // speed of light in vacuum [mm/ps]
    BFT refind = 1.4; // refractive index
    BFT alpha = A_Keijzer(refind); // boundary term
    BFT c = c0/refind;       // speed of light in medium [mm/ps]
    BFT freq = 100*1e6; // modulation frequency [Hz]
    BFT omega = 2.0*M_PI * freq*1e-12; // modulation frequency [cycles/ps]

    // Physical parameters, outer region
    BFT mua1 = 0.01; // absorption coefficient
    BFT mus1 = 1.0;  // scattering coefficient
    BFT kappa1 = 1.0/(3.0*(mua1+mus1));   // diffusion coefficient
    RT waveNumber1 = sqrt (RT(mua1/kappa1, omega/(c*kappa1))); // outer region

    // Physical parameters, inner region
    BFT mua2 = 0.02; // absorption coefficient
    BFT mus2 = 0.5;  // scattering coefficient
    BFT kappa2 = 1.0/(3.0*(mua2+mus2));   // diffusion coefficient
    RT waveNumber2 = sqrt (RT(mua2/kappa2, omega/(c*kappa2))); // outer region

    // Physical parameters, inner region
    BFT mua3 = 0.005; // absorption coefficient
    BFT mus3 = 2.0;  // scattering coefficient
    BFT kappa3 = 1.0/(3.0*(mua3+mus3));   // diffusion coefficient
    RT waveNumber3 = sqrt (RT(mua3/kappa3, omega/(c*kappa3))); // outer region

    // Create sphere meshes on the fly
    std::auto_ptr<Grid> grid1 = CreateSphere(25.0, 1.0);
    std::auto_ptr<Grid> grid2 = CreateSphere(15.0, 1.0);
    std::auto_ptr<Grid> grid3 = CreateSphere(10.0, 1.0);

    // Initialize the spaces

    PiecewiseLinearContinuousScalarSpace<BFT> HplusHalfSpace1(*grid1);
    PiecewiseLinearContinuousScalarSpace<BFT> HplusHalfSpace2(*grid2);
    PiecewiseLinearContinuousScalarSpace<BFT> HplusHalfSpace3(*grid3);

    HplusHalfSpace1.assignDofs();
    HplusHalfSpace2.assignDofs();
    HplusHalfSpace3.assignDofs();

    // Define some default options.

    AssemblyOptions assemblyOptions;

    // We want to use ACA

    AcaOptions acaOptions; // Default parameters for ACA
    acaOptions.eps = 1e-5;
    assemblyOptions.switchToAcaMode(acaOptions);

    // Define the standard integration factory

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(1);
    NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);
    Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

    // We need the single layer, double layer, and the identity operator

    // mesh1 x mesh1
    BoundaryOperator<BFT, RT> slp11 =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace1),
                SHARED(HplusHalfSpace1), SHARED(HplusHalfSpace1), waveNumber1);
    BoundaryOperator<BFT, RT> dlp11 =
            modifiedHelmholtz3dDoubleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace1),
                SHARED(HplusHalfSpace1), SHARED(HplusHalfSpace1), waveNumber1);
    BoundaryOperator<BFT, RT> id11 =
            identityOperator<BFT, RT>(
                SHARED(context), SHARED(HplusHalfSpace1),
                SHARED(HplusHalfSpace1), SHARED(HplusHalfSpace1));

    // mesh2 x mesh2, wavenumber 1
    BoundaryOperator<BFT, RT> slp22_w1 =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace2),
                SHARED(HplusHalfSpace2), SHARED(HplusHalfSpace2), waveNumber1);
    BoundaryOperator<BFT, RT> dlp22_w1 =
            modifiedHelmholtz3dDoubleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace2),
                SHARED(HplusHalfSpace2), SHARED(HplusHalfSpace2), waveNumber1);
    BoundaryOperator<BFT, RT> id22 =
            identityOperator<BFT, RT>(
                SHARED(context), SHARED(HplusHalfSpace2),
                SHARED(HplusHalfSpace2), SHARED(HplusHalfSpace2));

    // mesh2 x mesh2, wavenumber 2
    BoundaryOperator<BFT, RT> slp22_w2 =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace2),
                SHARED(HplusHalfSpace2), SHARED(HplusHalfSpace2), waveNumber2);
    BoundaryOperator<BFT, RT> dlp22_w2 =
            modifiedHelmholtz3dDoubleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace2),
                SHARED(HplusHalfSpace2), SHARED(HplusHalfSpace2), waveNumber2);

    // mesh1 x mesh2
    BoundaryOperator<BFT, RT> slp12 =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace2),
                SHARED(HplusHalfSpace1), SHARED(HplusHalfSpace1), waveNumber1);
    BoundaryOperator<BFT, RT> dlp12 =
            modifiedHelmholtz3dDoubleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace2),
                SHARED(HplusHalfSpace1), SHARED(HplusHalfSpace1), waveNumber1);

    // mesh2 x mesh1
    BoundaryOperator<BFT, RT> slp21 =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace1),
                SHARED(HplusHalfSpace2), SHARED(HplusHalfSpace2), waveNumber1);
    BoundaryOperator<BFT, RT> dlp21 =
            modifiedHelmholtz3dDoubleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace1),
                SHARED(HplusHalfSpace2), SHARED(HplusHalfSpace2), waveNumber1);

    // mesh2 x mesh3
    BoundaryOperator<BFT, RT> slp23 =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace3),
                SHARED(HplusHalfSpace2), SHARED(HplusHalfSpace2), waveNumber2);
    BoundaryOperator<BFT, RT> dlp23 =
            modifiedHelmholtz3dDoubleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace3),
                SHARED(HplusHalfSpace2), SHARED(HplusHalfSpace2), waveNumber2);

    // mesh3 x mesh2
    BoundaryOperator<BFT, RT> slp32 =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace2),
                SHARED(HplusHalfSpace3), SHARED(HplusHalfSpace3), waveNumber2);
    BoundaryOperator<BFT, RT> dlp32 =
            modifiedHelmholtz3dDoubleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace2),
                SHARED(HplusHalfSpace3), SHARED(HplusHalfSpace3), waveNumber2);

    // mesh3 x mesh3, wavenumber 2
    BoundaryOperator<BFT, RT> slp33_w2 =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace3),
                SHARED(HplusHalfSpace3), SHARED(HplusHalfSpace3), waveNumber2);
    BoundaryOperator<BFT, RT> dlp33_w2 =
            modifiedHelmholtz3dDoubleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace3),
                SHARED(HplusHalfSpace3), SHARED(HplusHalfSpace3), waveNumber2);

    // mesh3 x mesh3, wavenumber 3
    BoundaryOperator<BFT, RT> slp33_w3 =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace3),
                SHARED(HplusHalfSpace3), SHARED(HplusHalfSpace3), waveNumber3);
    BoundaryOperator<BFT, RT> dlp33_w3 =
            modifiedHelmholtz3dDoubleLayerBoundaryOperator<BFT, RT, RT>(
                SHARED(context), SHARED(HplusHalfSpace3),
                SHARED(HplusHalfSpace3), SHARED(HplusHalfSpace3), waveNumber3);
    BoundaryOperator<BFT, RT> id33 =
            identityOperator<BFT, RT>(
                SHARED(context), SHARED(HplusHalfSpace3),
                SHARED(HplusHalfSpace3), SHARED(HplusHalfSpace3));

    BFT scale = 1.0/(2.0*alpha*kappa1);
    BoundaryOperator<BFT, RT> lhs_k11 = 0.5*id11 + dlp11 + scale*slp11;
    BoundaryOperator<BFT, RT> lhs_k12 = -1.0*dlp12; // sign flipped to accommodate normal direction
    BoundaryOperator<BFT, RT> lhs_k13 = -(1.0/kappa1)*slp12;
    BoundaryOperator<BFT, RT> lhs_k21 = dlp21 + scale*slp21;
    BoundaryOperator<BFT, RT> lhs_k22 = 0.5*id22 - dlp22_w1; // sign flipped to accommodate normal direction
    BoundaryOperator<BFT, RT> lhs_k23 = -(1.0/kappa1)*slp22_w1;
    BoundaryOperator<BFT, RT> lhs_k32 = 0.5*id22 + dlp22_w2;
    BoundaryOperator<BFT, RT> lhs_k33 = (1.0/kappa2)*slp22_w2;
    BoundaryOperator<BFT, RT> lhs_k34 = -1.0*dlp23;
    BoundaryOperator<BFT, RT> lhs_k35 = -(1.0/kappa2)*slp23;
    BoundaryOperator<BFT, RT> lhs_k42 = dlp32;
    BoundaryOperator<BFT, RT> lhs_k43 = (1.0/kappa2)*slp32;
    BoundaryOperator<BFT, RT> lhs_k44 = 0.5*id33 - dlp33_w2;
    BoundaryOperator<BFT, RT> lhs_k45 = -(1.0/kappa2)*slp33_w2;
    BoundaryOperator<BFT, RT> lhs_k54 = 0.5*id33 + dlp33_w3;
    BoundaryOperator<BFT, RT> lhs_k55 = (1.0/kappa3)*slp33_w3;

    BlockedOperatorStructure<BFT, RT> structure;
    structure.setBlock(0, 0, lhs_k11);
    structure.setBlock(0, 1, lhs_k12);
    structure.setBlock(0, 2, lhs_k13);
    structure.setBlock(1, 0, lhs_k21);
    structure.setBlock(1, 1, lhs_k22);
    structure.setBlock(1, 2, lhs_k23);
    structure.setBlock(2, 1, lhs_k32);
    structure.setBlock(2, 2, lhs_k33);
    structure.setBlock(2, 3, lhs_k34);
    structure.setBlock(2, 4, lhs_k35);
    structure.setBlock(3, 1, lhs_k42);
    structure.setBlock(3, 2, lhs_k43);
    structure.setBlock(3, 3, lhs_k44);
    structure.setBlock(3, 4, lhs_k45);
    structure.setBlock(4, 3, lhs_k54);
    structure.setBlock(4, 4, lhs_k55);
    BlockedBoundaryOperator<BFT, RT> blockedOp(structure);

    // Grid functions for the RHS

    // TODO: remove the necessity of creating "dummy" functions
    // corresponding to zero blocks.
    BoundaryOperator<BFT, RT> rhs1 = scale*slp11;
    BoundaryOperator<BFT, RT> rhs2 = scale*slp21;

    std::vector<GridFunction<BFT, RT> > blockedRhs(5);
    blockedRhs[0] = rhs1 * GridFunction<BFT, RT>(
                SHARED(context),
                SHARED(HplusHalfSpace1), SHARED(HplusHalfSpace1),
                surfaceNormalIndependentFunction(MyFunctor(waveNumber1)));
    blockedRhs[1] = rhs2 * GridFunction<BFT, RT>(
                SHARED(context),
                SHARED(HplusHalfSpace1), SHARED(HplusHalfSpace1),
                surfaceNormalIndependentFunction(MyFunctor(waveNumber1)));
    blockedRhs[2] = GridFunction<BFT, RT>(
                SHARED(context),
                SHARED(HplusHalfSpace2), SHARED(HplusHalfSpace2),
                surfaceNormalIndependentFunction(NullFunctor()));
    blockedRhs[3] = GridFunction<BFT, RT>(
                SHARED(context),
                SHARED(HplusHalfSpace3), SHARED(HplusHalfSpace3),
                surfaceNormalIndependentFunction(NullFunctor()));
    blockedRhs[4] = GridFunction<BFT, RT>(
                SHARED(context),
                SHARED(HplusHalfSpace3), SHARED(HplusHalfSpace3),
                surfaceNormalIndependentFunction(NullFunctor()));

    // Initialize the solver

    const double solverTol = 1e-10;
    DefaultIterativeSolver<BFT, RT> solver(blockedOp);
    solver.initializeSolver(defaultGmresParameterList(solverTol));

    // Solve

    BlockedSolution<BFT, RT> solution = solver.solve(blockedRhs);

    std::cout << solution.solverMessage() << std::endl;
    arma::Col<RT> solutionVectorBlock1 = solution.gridFunction(0).coefficients();
    arma::Col<RT> solutionVectorBlock2 = solution.gridFunction(1).coefficients();
    arma::Col<RT> solutionVectorBlock3 = solution.gridFunction(2).coefficients();
    arma::Col<RT> solutionVectorBlock4 = solution.gridFunction(3).coefficients();
    arma::Col<RT> solutionVectorBlock5 = solution.gridFunction(4).coefficients();

    arma::diskio::save_raw_ascii(solutionVectorBlock1, "sol1.txt");
    arma::diskio::save_raw_ascii(solutionVectorBlock2, "sol2.txt");
    arma::diskio::save_raw_ascii(solutionVectorBlock3, "sol3.txt");
    arma::diskio::save_raw_ascii(solutionVectorBlock4, "sol4.txt");
    arma::diskio::save_raw_ascii(solutionVectorBlock5, "sol5.txt");

    solution.gridFunction(0).exportToVtk(VtkWriter::VERTEX_DATA, "gf1", "gf1");
    solution.gridFunction(1).exportToVtk(VtkWriter::VERTEX_DATA, "gf2", "gf2");
    solution.gridFunction(2).exportToVtk(VtkWriter::VERTEX_DATA, "gf3", "gf3");
    solution.gridFunction(3).exportToVtk(VtkWriter::VERTEX_DATA, "gf4", "gf4");
    solution.gridFunction(4).exportToVtk(VtkWriter::VERTEX_DATA, "gf5", "gf5");
}
