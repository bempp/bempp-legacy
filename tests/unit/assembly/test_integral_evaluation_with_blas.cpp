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

#include "bempp/common/config_ahmed.hpp"

#ifdef WITH_AHMED

#include "../type_template.hpp"
#include "../type_combinations.hpp"
#include "../check_arrays_are_close.hpp"

#include "assembly/context.hpp"
#include "assembly/discrete_boundary_operator.hpp"
#include "assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_hypersingular_boundary_operator.hpp"
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/helmholtz_3d_hypersingular_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_double_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_hypersingular_boundary_operator.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"

#include "grid/grid_factory.hpp"

#include "space/piecewise_constant_scalar_space.hpp"
#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_linear_discontinuous_scalar_space.hpp"
#include "space/piecewise_polynomial_continuous_scalar_space.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/type_traits/is_complex.hpp>

using namespace Bempp;

namespace Bempp
{

namespace
{

template <typename T> T initWaveNumber();
template <> float initWaveNumber() { return 1.2f; }
template <> double initWaveNumber(){ return 1.2; }
template <> std::complex<float> initWaveNumber()
{ return std::complex<float>(1.2f, 0.7f); }
template <> std::complex<double> initWaveNumber()
{ return std::complex<double>(1.2, 0.7); }

}

}

// Tests

BOOST_AUTO_TEST_SUITE(IntegralEvaluationWithBlas)

BOOST_AUTO_TEST_CASE_TEMPLATE(blas_works_for_modified_helmholtz_3d_double_layer_operator_spaces_a,
                              Traits, basis_kernel_result_combinations)
{
    typedef typename Traits::BasisFunctionType BFT;
    typedef typename Traits::KernelType KT;
    typedef typename Traits::ResultType RT;
    typedef typename ScalarTraits<RT>::RealType RealType;

    KT waveNumber = initWaveNumber<KT>();

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
//        params, "../../meshes/sphere-h-0.4.msh", false /* verbose */);
    params, "meshes/cube-12-reoriented.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseQuads(
        new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 2));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsNoBlas;
    assemblyOptionsNoBlas.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextNoBlas(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsNoBlas));

    BoundaryOperator<BFT, RT> opNoBlas =
            modifiedHelmholtz3dDoubleLayerBoundaryOperator<BFT, KT, RT>(
                contextNoBlas, pwiseQuads, pwiseQuads, pwiseLinears, waveNumber);
    arma::Mat<RT> weakFormNoBlas = opNoBlas.weakForm()->asMatrix();

    AssemblyOptions assemblyOptionsBlas;
    assemblyOptionsBlas.enableBlasInQuadrature();
    assemblyOptionsBlas.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextBlas(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsBlas));

    BoundaryOperator<BFT, RT> opBlas =
            modifiedHelmholtz3dDoubleLayerBoundaryOperator<BFT, KT, RT>(
                contextBlas, pwiseQuads, pwiseQuads, pwiseLinears, waveNumber);
    arma::Mat<RT> weakFormBlas = opBlas.weakForm()->asMatrix();

    BOOST_CHECK(check_arrays_are_close<RT>(
                    weakFormNoBlas, weakFormBlas,
                    100 * std::numeric_limits<RealType>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(blas_works_for_modified_helmholtz_3d_double_layer_operator_spaces_b,
                              Traits, basis_kernel_result_combinations)
{
    typedef typename Traits::BasisFunctionType BFT;
    typedef typename Traits::KernelType KT;
    typedef typename Traits::ResultType RT;
    typedef typename ScalarTraits<RT>::RealType RealType;

    KT waveNumber = initWaveNumber<KT>();

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
//        params, "../../meshes/sphere-h-0.4.msh", false /* verbose */);
    params, "meshes/cube-12-reoriented.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseQuads(
        new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 2));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsNoBlas;
    assemblyOptionsNoBlas.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextNoBlas(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsNoBlas));

    BoundaryOperator<BFT, RT> opNoBlas =
            modifiedHelmholtz3dDoubleLayerBoundaryOperator<BFT, KT, RT>(
                contextNoBlas, pwiseLinears, pwiseLinears, pwiseQuads, waveNumber);
    arma::Mat<RT> weakFormNoBlas = opNoBlas.weakForm()->asMatrix();

    AssemblyOptions assemblyOptionsBlas;
    assemblyOptionsBlas.enableBlasInQuadrature();
    assemblyOptionsBlas.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextBlas(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsBlas));

    BoundaryOperator<BFT, RT> opBlas =
            modifiedHelmholtz3dDoubleLayerBoundaryOperator<BFT, KT, RT>(
                contextBlas, pwiseLinears, pwiseLinears, pwiseQuads, waveNumber);
    arma::Mat<RT> weakFormBlas = opBlas.weakForm()->asMatrix();

    BOOST_CHECK(check_arrays_are_close<RT>(
                    weakFormNoBlas, weakFormBlas,
                    100 * std::numeric_limits<RealType>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(blas_works_for_laplace_3d_hypersingular_operator_spaces_a,
                              Traits, basis_result_combinations)
{
    typedef typename Traits::BasisFunctionType BFT;
    typedef typename Traits::ResultType RT;
    typedef typename ScalarTraits<RT>::RealType RealType;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
    //        params, "../../meshes/sphere-h-0.4.msh", false /* verbose */);
                params, "meshes/cube-12-reoriented.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseQuads(
        new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 2));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsNoBlas;
    assemblyOptionsNoBlas.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextNoBlas(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsNoBlas));

    BoundaryOperator<BFT, RT> opNoBlas =
            laplace3dHypersingularBoundaryOperator<BFT, RT>(
                contextNoBlas, pwiseQuads, pwiseQuads, pwiseLinears);
    arma::Mat<RT> weakFormNoBlas = opNoBlas.weakForm()->asMatrix();

    AssemblyOptions assemblyOptionsBlas;
    assemblyOptionsBlas.enableBlasInQuadrature();
    assemblyOptionsBlas.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextBlas(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsBlas));

    BoundaryOperator<BFT, RT> opBlas =
            laplace3dHypersingularBoundaryOperator<BFT, RT>(
                contextBlas, pwiseQuads, pwiseQuads, pwiseLinears);
    arma::Mat<RT> weakFormBlas = opBlas.weakForm()->asMatrix();

    BOOST_CHECK(check_arrays_are_close<RT>(
                    weakFormNoBlas, weakFormBlas,
                    100 * std::numeric_limits<RealType>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(blas_works_for_laplace_3d_hypersingular_operator_spaces_b,
                              Traits, basis_result_combinations)
{
    typedef typename Traits::BasisFunctionType BFT;
    typedef typename Traits::ResultType RT;
    typedef typename ScalarTraits<RT>::RealType RealType;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
    //        params, "../../meshes/sphere-h-0.4.msh", false /* verbose */);
                params, "meshes/cube-12-reoriented.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseQuads(
        new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 2));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsNoBlas;
    assemblyOptionsNoBlas.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextNoBlas(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsNoBlas));

    BoundaryOperator<BFT, RT> opNoBlas =
            laplace3dHypersingularBoundaryOperator<BFT, RT>(
                contextNoBlas, pwiseLinears, pwiseQuads, pwiseQuads);
    arma::Mat<RT> weakFormNoBlas = opNoBlas.weakForm()->asMatrix();

    AssemblyOptions assemblyOptionsBlas;
    assemblyOptionsBlas.enableBlasInQuadrature();
    assemblyOptionsBlas.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextBlas(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsBlas));

    BoundaryOperator<BFT, RT> opBlas =
            laplace3dHypersingularBoundaryOperator<BFT, RT>(
                contextBlas, pwiseLinears, pwiseQuads, pwiseQuads);
    arma::Mat<RT> weakFormBlas = opBlas.weakForm()->asMatrix();

    BOOST_CHECK(check_arrays_are_close<RT>(
                    weakFormNoBlas, weakFormBlas,
                    100 * std::numeric_limits<RealType>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(blas_works_for_modified_helmholtz_3d_hypersingular_operator_spaces_a,
                              Traits, basis_kernel_result_combinations)
{
    typedef typename Traits::BasisFunctionType BFT;
    typedef typename Traits::KernelType KT;
    typedef typename Traits::ResultType RT;
    typedef typename ScalarTraits<RT>::RealType RealType;

    KT waveNumber = initWaveNumber<KT>();

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
//        params, "../../meshes/sphere-h-0.4.msh", false /* verbose */);
    params, "meshes/cube-12-reoriented.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseQuads(
        new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 2));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsNoBlas;
    assemblyOptionsNoBlas.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextNoBlas(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsNoBlas));

    BoundaryOperator<BFT, RT> opNoBlas =
            modifiedHelmholtz3dHypersingularBoundaryOperator<BFT, KT, RT>(
                contextNoBlas, pwiseQuads, pwiseQuads, pwiseLinears, waveNumber);
    arma::Mat<RT> weakFormNoBlas = opNoBlas.weakForm()->asMatrix();

    AssemblyOptions assemblyOptionsBlas;
    assemblyOptionsBlas.enableBlasInQuadrature();
    assemblyOptionsBlas.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextBlas(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsBlas));

    BoundaryOperator<BFT, RT> opBlas =
            modifiedHelmholtz3dHypersingularBoundaryOperator<BFT, KT, RT>(
                contextBlas, pwiseQuads, pwiseQuads, pwiseLinears, waveNumber);
    arma::Mat<RT> weakFormBlas = opBlas.weakForm()->asMatrix();

    BOOST_CHECK(check_arrays_are_close<RT>(
                    weakFormNoBlas, weakFormBlas,
                    100 * std::numeric_limits<RealType>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(blas_works_for_modified_helmholtz_3d_hypersingular_operator_spaces_b,
                              Traits, basis_kernel_result_combinations)
{
    typedef typename Traits::BasisFunctionType BFT;
    typedef typename Traits::KernelType KT;
    typedef typename Traits::ResultType RT;
    typedef typename ScalarTraits<RT>::RealType RealType;

    KT waveNumber = initWaveNumber<KT>();

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
//        params, "../../meshes/sphere-h-0.4.msh", false /* verbose */);
    params, "meshes/cube-12-reoriented.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseQuads(
        new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 2));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsNoBlas;
    assemblyOptionsNoBlas.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextNoBlas(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsNoBlas));

    BoundaryOperator<BFT, RT> opNoBlas =
            modifiedHelmholtz3dHypersingularBoundaryOperator<BFT, KT, RT>(
                contextNoBlas, pwiseLinears, pwiseLinears, pwiseQuads, waveNumber);
    arma::Mat<RT> weakFormNoBlas = opNoBlas.weakForm()->asMatrix();

    AssemblyOptions assemblyOptionsBlas;
    assemblyOptionsBlas.enableBlasInQuadrature();
    assemblyOptionsBlas.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextBlas(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsBlas));

    BoundaryOperator<BFT, RT> opBlas =
            modifiedHelmholtz3dHypersingularBoundaryOperator<BFT, KT, RT>(
                contextBlas, pwiseLinears, pwiseLinears, pwiseQuads, waveNumber);
    arma::Mat<RT> weakFormBlas = opBlas.weakForm()->asMatrix();

    BOOST_CHECK(check_arrays_are_close<RT>(
                    weakFormNoBlas, weakFormBlas,
                    100 * std::numeric_limits<RealType>::epsilon()));
}

BOOST_AUTO_TEST_SUITE_END()

#endif // WITH_AHMED
