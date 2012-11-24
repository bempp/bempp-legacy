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

#include "../check_arrays_are_close.hpp"
#include "../type_template.hpp"

#include "assembly/assembly_options.hpp"
#include "assembly/discrete_boundary_operator.hpp"
#include "assembly/boundary_operator.hpp"
#include "assembly/context.hpp"
#include "assembly/identity_operator.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"

#include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"

#include "common/boost_make_shared_fwd.hpp"

#include "grid/grid_factory.hpp"
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

template <typename T>
T initMultiplier();

template <>
float initMultiplier()
{
    return 1.2f;
}

template <>
double initMultiplier()
{
    return 1.2;
}

template <>
std::complex<float> initMultiplier()
{
    return std::complex<float>(1.2f, 0.7f);
}

template <>
std::complex<double> initMultiplier()
{
    return std::complex<double>(1.2, 0.7);
}

} // namespace

BOOST_AUTO_TEST_SUITE(ScaledDiscreteDenseBoundaryOperator)

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_alpha_equal_to_2_and_beta_equal_to_0_and_y_initialized_to_nans, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;

    const int dimGrid = 2;
    typedef double ctype;
    arma::Col<double> lowerLeft(dimGrid);
    arma::Col<double> upperRight(dimGrid);
    arma::Col<unsigned int> nElements(dimGrid);
    lowerLeft.fill(0);
    upperRight.fill(1);
    nElements(0) = 2;
    nElements(1) = 3;

    shared_ptr<Grid> grid =
        Bempp::GridFactory::createStructuredGrid(
            params, lowerLeft, upperRight, nElements);

    PiecewiseLinearContinuousScalarSpace<BFT> pwiseLinears(grid);
    PiecewiseConstantScalarSpace<BFT> pwiseConstants(grid);

    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    NumericalQuadratureStrategy<BFT, RT> quadStrategy;

    Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

    const RT multiplier = initMultiplier<RT>();

    BoundaryOperator<BFT, RT> id =
        identityOperator<BFT, RT>(
        make_shared_from_ref(context),
        make_shared_from_ref(pwiseConstants),
        make_shared_from_ref(pwiseConstants),
        make_shared_from_ref(pwiseLinears));

    BoundaryOperator<BFT, RT> op = multiplier * id;

    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = op.weakForm();

    RT alpha = static_cast<RT>(2.);
    RT beta = static_cast<RT>(0.);

    arma::Col<RT> x(dop->columnCount());
    x.fill(1.);
    arma::Col<RT> y(dop->rowCount());
    y.fill(std::numeric_limits<CT>::quiet_NaN());

    arma::Col<RT> expected = alpha * dop->asMatrix() * x;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(y.is_finite());
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_SUITE_END()
