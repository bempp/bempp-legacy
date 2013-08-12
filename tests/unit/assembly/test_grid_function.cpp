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
#include "assembly/boundary_operator.hpp"
#include "assembly/context.hpp"
#include "assembly/grid_function.hpp"
#include "assembly/identity_operator.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "assembly/surface_normal_independent_function.hpp"

#include "common/scalar_traits.hpp"

#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"

#include <boost/test/floating_point_comparison.hpp>

using namespace Bempp;

template <typename ValueType_>
class ConstantFunction
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    int argumentDimension() const { return 3; }
    int resultDimension() const { return 1; }

    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        result(0) = 2.;
    }
};

template <typename ValueType_>
class SinusoidalFunction
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    int argumentDimension() const { return 3; }
    int resultDimension() const { return 1; }

    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        result(0) = point(0);
    }
};

template <typename ValueType_>
class ExponentialFunction
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    int argumentDimension() const { return 3; }
    int resultDimension() const { return 1; }

    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        result(0) = std::exp(ValueType(0., point(0)));
    }
};

// Tests

BOOST_AUTO_TEST_SUITE(GridFunction)

BOOST_AUTO_TEST_CASE_TEMPLATE(L2Norm_works_for_constant_function_and_piecewise_constants, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../examples/meshes/sphere-h-0.1.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        new PiecewiseConstantScalarSpace<BFT>(grid));

    AccuracyOptions accuracyOptions;
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));
    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    Bempp::GridFunction<BFT, RT> fun(context, space, space,
                surfaceNormalIndependentFunction(
                    ConstantFunction<RT>()));

    CT norm = fun.L2Norm();
    CT expectedNorm = sqrt(2. * 2. * 4. * M_PI);
    BOOST_CHECK_CLOSE(norm, expectedNorm, 1 /* percent */);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(L2Norm_works_for_sinusoidal_function_and_piecewise_constants, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../examples/meshes/sphere-h-0.1.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        new PiecewiseConstantScalarSpace<BFT>(grid));

    AccuracyOptions accuracyOptions;
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));
    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    Bempp::GridFunction<BFT, RT> fun(context, space, space,
                surfaceNormalIndependentFunction(
                    SinusoidalFunction<RT>()));
    exportToVtk(fun, VtkWriter::CELL_DATA, "fun", "fun");

    CT norm = fun.L2Norm();
    CT expectedNorm = sqrt(4. / 3. * M_PI);
    BOOST_CHECK_CLOSE(norm, expectedNorm, 1 /* percent */);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(L2Norm_works_for_exponential_function_and_piecewise_constants, ResultType, complex_result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../examples/meshes/sphere-h-0.1.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        new PiecewiseConstantScalarSpace<BFT>(grid));

    AccuracyOptions accuracyOptions;
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));
    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    Bempp::GridFunction<BFT, RT> fun(context, space, space,
                surfaceNormalIndependentFunction(
                    ExponentialFunction<RT>()));

    CT norm = fun.L2Norm();
    CT expectedNorm = sqrt(4 * M_PI);
    BOOST_CHECK_CLOSE(norm, expectedNorm, 1 /* percent */);
}

BOOST_AUTO_TEST_SUITE_END()
