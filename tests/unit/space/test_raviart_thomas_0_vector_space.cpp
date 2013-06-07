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

#include "../assembly/create_regular_grid.hpp"

#include "../check_arrays_are_close.hpp"
#include "../type_template.hpp"

#include "assembly/assembly_options.hpp"
#include "assembly/context.hpp"
#include "assembly/grid_function.hpp"
#include "assembly/l2_norm.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "assembly/surface_normal_independent_function.hpp"

#include "common/scalar_traits.hpp"

#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"

#include "space/raviart_thomas_0_vector_space.hpp"

#include <boost/type_traits/is_complex.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace Bempp;

/* Grid:

  -----------------
  |  /|  /|  /|  /|
  | / | / | / | / |
  |/  |/  |/  |/  |
  -----------------
  |  /|  /|  /|  /|
  | / | / | / | / |
  |/  |/  |/  |/  |
  -----------------
*/

// This function is composed of all the basis functions having a nonzero
// normal to the diagonal edges of the grid.

template <typename ValueType_>
class MyFunction
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    int argumentDimension() const { return 3; }
    int resultDimension() const { return 3; }

    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        CoordinateType x = point(0), y = point(1);
        CoordinateType red_x = fmod(x, 0.2);
        CoordinateType red_y = fmod(y, 0.2);
        if (red_x  > red_y) {
            result(0) = 0.2 - red_x;
            result(1) = -red_y;
        } else {
            result(0) = red_x;
            result(1) = -0.2 + red_y;
        }
        result(2) = 0.;
    }
};

// Tests

BOOST_AUTO_TEST_SUITE(RaviartThomas0VectorSpace_)

BOOST_AUTO_TEST_CASE_TEMPLATE(my_function_can_be_expanded_in_rt_space_with_dofs_on_boundary, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    shared_ptr<Grid> grid = createRegularTriangularGrid(5, 10, 1., 2.);

    shared_ptr<Space<BFT> > space(
        new RaviartThomas0VectorSpace<BFT>(grid,
                                           true /* put DOFs on boundary */));

    AccuracyOptions accuracyOptions;
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
        new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));
    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    GridFunction<BFT, RT> function(
        context, space, space,
        surfaceNormalIndependentFunction(
                    MyFunction<RT>()));
    CT absoluteError, relativeError;
    estimateL2Error(
        function, surfaceNormalIndependentFunction(
                    MyFunction<RT>()),
        *quadStrategy, absoluteError, relativeError);
    BOOST_CHECK_SMALL(relativeError,
                      1000 * std::numeric_limits<CT>::epsilon() /* percent */);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(my_function_can_be_expanded_in_rt_space_without_dofs_on_boundary, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    shared_ptr<Grid> grid = createRegularTriangularGrid(5, 10, 1., 2.);

    shared_ptr<Space<BFT> > space(
        new RaviartThomas0VectorSpace<BFT>(grid,
                                           false /* put DOFs on boundary */));

    AccuracyOptions accuracyOptions;
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
        new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));
    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    GridFunction<BFT, RT> function(
        context, space, space,
        surfaceNormalIndependentFunction(
                    MyFunction<RT>()));
    CT absoluteError, relativeError;
    estimateL2Error(
        function, surfaceNormalIndependentFunction(
                    MyFunction<RT>()),
        *quadStrategy, absoluteError, relativeError);
    BOOST_CHECK_SMALL(relativeError,
                      1000 * std::numeric_limits<CT>::epsilon() /* percent */);
}

BOOST_AUTO_TEST_SUITE_END()
