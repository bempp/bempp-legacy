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

#include "../check_arrays_are_close.hpp"
#include "../type_template.hpp"
#include "../random_arrays.hpp"

#include "create_regular_grid.hpp"

#include "assembly/assembly_options.hpp"
#include "assembly/discrete_aca_boundary_operator.hpp"
#include "assembly/discrete_boundary_operator.hpp"
#include "assembly/boundary_operator.hpp"
#include "assembly/context.hpp"
#include "assembly/identity_operator.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"

#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_hypersingular_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"

#include "grid/grid.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"

#include <algorithm>
#include "common/armadillo_fwd.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/version.hpp>
#include <complex>

// Tests

using namespace Bempp;

namespace 
{

template <typename T> T initWaveNumber();
template <> float initWaveNumber() { return 1.2f; }
template <> double initWaveNumber(){ return 1.2; }
template <> std::complex<float> initWaveNumber() 
{ return std::complex<float>(1.2f, 0.7f); }
template <> std::complex<double> initWaveNumber()
{ return std::complex<double>(1.2, 0.7); }

template <typename BFT, typename RT> 
struct DiscreteAcaBoundaryOperatorFixture
{
    DiscreteAcaBoundaryOperatorFixture()
    {
        grid = createRegularTriangularGrid(4, 7);

        shared_ptr<Space<BFT> > pwiseConstants(
            new PiecewiseConstantScalarSpace<BFT>(grid));
        shared_ptr<Space<BFT> > pwiseLinears(
            new PiecewiseLinearContinuousScalarSpace<BFT>(grid));

        AssemblyOptions assemblyOptions;
        assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
        AcaOptions acaOptions;
        acaOptions.minimumBlockSize = 2;
        assemblyOptions.switchToAcaMode(acaOptions);
        shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy( 
            new NumericalQuadratureStrategy<BFT, RT>);
        shared_ptr<Context<BFT, RT> > context(
            new Context<BFT, RT>(quadStrategy, assemblyOptions));
        
        const RT waveNumber = initWaveNumber<RT>();
        
        op = modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, RT, RT>(
            context, pwiseConstants, pwiseConstants, pwiseLinears, waveNumber);
    }

    shared_ptr<Grid> grid;
    BoundaryOperator<BFT, RT> op;
};

template <typename BFT, typename RT>
struct DiscreteHermitianAcaBoundaryOperatorFixture
{
    DiscreteHermitianAcaBoundaryOperatorFixture()
    {
        grid = createRegularTriangularGrid(4, 7);

        shared_ptr<Space<BFT> > pwiseConstants(
            new PiecewiseConstantScalarSpace<BFT>(grid));

        AssemblyOptions assemblyOptions;
        assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
        AcaOptions acaOptions;
        acaOptions.minimumBlockSize = 2;
        assemblyOptions.switchToAcaMode(acaOptions);
        shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
            new NumericalQuadratureStrategy<BFT, RT>);
        shared_ptr<Context<BFT, RT> > context(
            new Context<BFT, RT>(quadStrategy, assemblyOptions));

        op = laplace3dSingleLayerBoundaryOperator<BFT, RT>(
            context, pwiseConstants, pwiseConstants, pwiseConstants, "SLP", HERMITIAN);
    }

    shared_ptr<Grid> grid;
    BoundaryOperator<BFT, RT> op;
};

} // namespace

BOOST_AUTO_TEST_SUITE(DiscreteAcaBoundaryOperator)

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_alpha_equal_to_2_and_beta_equal_to_0_and_y_initialized_to_nans, ResultType, result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    DiscreteAcaBoundaryOperatorFixture<BFT, RT> fixture;
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.op.weakForm();

    RT alpha = static_cast<RT>(2.);
    RT beta = static_cast<RT>(0.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y(dop->rowCount());
    y.fill(std::numeric_limits<CT>::quiet_NaN());

    arma::Col<RT> expected = alpha * dop->asMatrix() * x;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);
    
    BOOST_CHECK(y.is_finite());
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected, 
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_alpha_equal_to_2_plus_3j_and_beta_equal_to_0_and_y_initialized_to_nans, ResultType, complex_result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    DiscreteAcaBoundaryOperatorFixture<BFT, RT> fixture;
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.op.weakForm();

    RT alpha = static_cast<RT>(2., 3.);
    RT beta = static_cast<RT>(0.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y(dop->rowCount());
    y.fill(std::numeric_limits<CT>::quiet_NaN());

    arma::Col<RT> expected = alpha * dop->asMatrix() * x;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);
    
    BOOST_CHECK(y.is_finite());
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected, 
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_alpha_equal_to_2_and_beta_equal_to_3, ResultType, result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    DiscreteAcaBoundaryOperatorFixture<BFT, RT> fixture;
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.op.weakForm();

    RT alpha = static_cast<RT>(2.);
    RT beta = static_cast<RT>(3.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y = generateRandomVector<RT>(dop->rowCount());

    arma::Col<RT> expected = alpha * dop->asMatrix() * x + beta * y;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);
    
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected, 
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_alpha_equal_to_2_plus_3j_and_beta_equal_to_4_minus_5j, ResultType, complex_result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    DiscreteAcaBoundaryOperatorFixture<BFT, RT> fixture;
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.op.weakForm();

    RT alpha = static_cast<RT>(2., 3.);
    RT beta = static_cast<RT>(4., -5.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y = generateRandomVector<RT>(dop->rowCount());

    arma::Col<RT> expected = alpha * dop->asMatrix() * x + beta * y;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);
    
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected, 
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_alpha_equal_to_2_and_beta_equal_to_3_and_conjugate_transpose, ResultType, result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    DiscreteAcaBoundaryOperatorFixture<BFT, RT> fixture;
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.op.weakForm();

    RT alpha = static_cast<RT>(2.);
    RT beta = static_cast<RT>(3.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->rowCount());
    arma::Col<RT> y = generateRandomVector<RT>(dop->columnCount());

    // .t() gives conjugate transpose for complex matrices
    arma::Col<RT> expected = alpha * dop->asMatrix().t() * x + beta * y;

    dop->apply(CONJUGATE_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_alpha_equal_to_2_plus_3j_and_beta_equal_to_4_minus_5j_and_conjugate_transpose, ResultType, complex_result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    DiscreteAcaBoundaryOperatorFixture<BFT, RT> fixture;
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.op.weakForm();

    RT alpha = static_cast<RT>(2., 3.);
    RT beta = static_cast<RT>(4., -5.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->rowCount());
    arma::Col<RT> y = generateRandomVector<RT>(dop->columnCount());

    // .t() gives conjugate transpose for complex matrices
    arma::Col<RT> expected = alpha * dop->asMatrix().t() * x + beta * y;

    dop->apply(CONJUGATE_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_alpha_equal_to_2_and_beta_equal_to_3_and_real_symmetric_operator, ResultType, result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    DiscreteHermitianAcaBoundaryOperatorFixture<BFT, RT> fixture;
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.op.weakForm();

    RT alpha = static_cast<RT>(2.);
    RT beta = static_cast<RT>(3.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->rowCount());
    arma::Col<RT> y = generateRandomVector<RT>(dop->columnCount());

    arma::Col<RT> expected = alpha * dop->asMatrix() * x + beta * y;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_alpha_equal_to_2_and_beta_equal_to_3_and_conjugate_transpose_and_real_symmetric_operator, ResultType, result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    DiscreteHermitianAcaBoundaryOperatorFixture<BFT, RT> fixture;
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.op.weakForm();

    RT alpha = static_cast<RT>(2.);
    RT beta = static_cast<RT>(3.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->rowCount());
    arma::Col<RT> y = generateRandomVector<RT>(dop->columnCount());

    // .t() gives conjugate transpose for complex matrices
    arma::Col<RT> expected = alpha * dop->asMatrix().t() * x + beta * y;

    dop->apply(CONJUGATE_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(acaOperatorSum_works_correctly_for_nonsymmetric_operators, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    shared_ptr<Grid> grid = createRegularTriangularGrid(10, 3);

    shared_ptr<Space<BFT> > pwiseConstants(
        new PiecewiseConstantScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));

    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    assemblyOptions.switchToAcaMode(AcaOptions());
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
        new NumericalQuadratureStrategy<BFT, RT>);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    const RT waveNumber = initWaveNumber<RT>();

    BoundaryOperator<BFT, RT> op =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, RT, RT>(
        context, pwiseConstants, pwiseConstants, pwiseLinears, waveNumber);
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = op.weakForm();

    BoundaryOperator<BFT, RT> op2 =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                context, pwiseConstants, pwiseConstants, pwiseLinears);
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop2 = op2.weakForm();

    arma::Mat<RT> expected = dop->asMatrix() + dop2->asMatrix();

    const double eps = 1e-4;
    shared_ptr<const DiscreteBoundaryOperator<RT> > acaSum =
            acaOperatorSum(dop,
                           dop2,
                           eps, UINT_MAX);

    BOOST_CHECK(check_arrays_are_close<RT>(acaSum->asMatrix(), expected,
                                           eps));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(acaOperatorSum_works_correctly_for_real_hermitian_operators,
                              ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    shared_ptr<Grid> grid = createRegularTriangularGrid(10, 3);

    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));

    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    assemblyOptions.switchToAcaMode(AcaOptions());
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
        new NumericalQuadratureStrategy<BFT, RT>);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    BoundaryOperator<BFT, RT> op =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                context, pwiseLinears, pwiseLinears, pwiseLinears, "SLP",
                SYMMETRIC | HERMITIAN);
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = op.weakForm();

    BoundaryOperator<BFT, RT> op2 =
            laplace3dHypersingularBoundaryOperator<BFT, RT>(
                context, pwiseLinears, pwiseLinears, pwiseLinears, "Hyp",
                SYMMETRIC | HERMITIAN);
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop2 = op2.weakForm();

    arma::Mat<RT> expected = dop->asMatrix() + dop2->asMatrix();

    const double eps = 1e-4;
    shared_ptr<const DiscreteBoundaryOperator<RT> > acaSum =
            acaOperatorSum(dop,
                           dop2,
                           eps, UINT_MAX);

    BOOST_CHECK(check_arrays_are_close<RT>(acaSum->asMatrix(), expected,
                                           eps));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(scaledAcaOperator_works_correctly_for_nonsymmetric_operators,
                              ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    shared_ptr<Grid> grid = createRegularTriangularGrid(10, 3);

    shared_ptr<Space<BFT> > pwiseConstants(
        new PiecewiseConstantScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));

    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    assemblyOptions.switchToAcaMode(AcaOptions());
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
        new NumericalQuadratureStrategy<BFT, RT>);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    const RT waveNumber = initWaveNumber<RT>();

    BoundaryOperator<BFT, RT> op =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, RT, RT>(
        context, pwiseConstants, pwiseConstants, pwiseLinears, waveNumber);
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = op.weakForm();

    arma::Mat<RT> expected = waveNumber * dop->asMatrix();

    shared_ptr<const DiscreteBoundaryOperator<RT> > scaled =
            scaledAcaOperator(dop,
                              waveNumber);

    BOOST_CHECK(check_arrays_are_close<RT>(scaled->asMatrix(), expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(scaledAcaOperator_works_correctly_for_real_hermitian_operators,
                              ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    shared_ptr<Grid> grid = createRegularTriangularGrid(10, 3);

    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));

    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    assemblyOptions.switchToAcaMode(AcaOptions());
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
        new NumericalQuadratureStrategy<BFT, RT>);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    BoundaryOperator<BFT, RT> op =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                context, pwiseLinears, pwiseLinears, pwiseLinears, "SLP",
                SYMMETRIC | HERMITIAN);
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = op.weakForm();

    const RT multiplier = 5.25;
    arma::Mat<RT> expected = multiplier * dop->asMatrix();

    shared_ptr<const DiscreteBoundaryOperator<RT> > scaled =
            scaledAcaOperator(dop,
                              multiplier);

    BOOST_CHECK(check_arrays_are_close<RT>(scaled->asMatrix(), expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_SUITE_END()

#endif // WITH_AHMED
