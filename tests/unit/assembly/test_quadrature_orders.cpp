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
#include "assembly/general_elementary_singular_integral_operator_imp.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"

#include "common/boost_make_shared_fwd.hpp"

#include "fiber/geometrical_data.hpp"
#include "fiber/scalar_function_value_functor.hpp"
#include "fiber/simple_test_scalar_kernel_trial_integrand_functor.hpp"

#include "grid/grid_factory.hpp"
#include "grid/grid.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_polynomial_continuous_scalar_space.hpp"
#include "space/piecewise_polynomial_discontinuous_scalar_space.hpp"

#include <algorithm>
#include "common/armadillo_fwd.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/version.hpp>
#include <complex>

// Tests

using namespace Bempp;

template <typename ValueType_>
class ConstantKernelFunctor
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    int kernelCount() const { return 1; }
    int kernelRowCount(int /* kernelIndex */) const { return 1; }
    int kernelColCount(int /* kernelIndex */) const { return 1; }

    void addGeometricalDependencies(size_t& testGeomDeps, size_t& trialGeomDeps) const {
    }

    template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
    void evaluate(
            const Fiber::ConstGeometricalDataSlice<CoordinateType>& testGeomData,
            const Fiber::ConstGeometricalDataSlice<CoordinateType>& trialGeomData,
            CollectionOf2dSlicesOfNdArrays<ValueType>& result) const {
        result[0](0, 0) = 1.;
    }
};

template <typename BasisFunctionType, typename KernelType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
constantBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY)
{
    typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

    typedef ConstantKernelFunctor<KernelType>
    KernelFunctor;
    typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
    TransformationFunctor;
    typedef Fiber::SimpleTestScalarKernelTrialIntegrandFunctorExt<
    BasisFunctionType, KernelType, ResultType, 1> IntegrandFunctor;

    typedef GeneralElementarySingularIntegralOperator<
            BasisFunctionType, KernelType, ResultType> Op;
    shared_ptr<Op> newOp(new Op(
                             domain, range, dualToRange, label, symmetry,
                             KernelFunctor(),
                             TransformationFunctor(),
                             TransformationFunctor(),
                             IntegrandFunctor()));
    return BoundaryOperator<BasisFunctionType, ResultType>(context, newOp);
}

BOOST_AUTO_TEST_SUITE(QuadratureOrders)

BOOST_AUTO_TEST_CASE_TEMPLATE(default_quadrature_order_is_sufficient_for_piecewise_linears,
                              BasisFunctionType, basis_function_types)
{
    typedef BasisFunctionType BFT;
    typedef BFT RT;
    typedef BFT KT;
    typedef typename Fiber::ScalarTraits<BFT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
                params, "meshes/sphere-ico-1.msh",
                false /* verbose */);

    PiecewiseLinearContinuousScalarSpace<BFT> space(grid);

    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    AccuracyOptions accuracyOptions;
    NumericalQuadratureStrategy<BFT, RT> defaultQuadStrategy(accuracyOptions);
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(5);
    accuracyOptions.doubleSingular.setRelativeQuadratureOrder(5);
    NumericalQuadratureStrategy<BFT, RT> accurateQuadStrategy(accuracyOptions);

    Context<BFT, RT> defaultContext(
                make_shared_from_ref(defaultQuadStrategy), assemblyOptions);
    Context<BFT, RT> accurateContext(
                make_shared_from_ref(accurateQuadStrategy), assemblyOptions);

    BoundaryOperator<BFT, RT> opDefault =
            constantBoundaryOperator<BFT, KT, RT>(
                make_shared_from_ref(defaultContext),
                make_shared_from_ref(space),
                make_shared_from_ref(space),
                make_shared_from_ref(space));
    BoundaryOperator<BFT, RT> opAccurate =
            constantBoundaryOperator<BFT, KT, RT>(
                make_shared_from_ref(accurateContext),
                make_shared_from_ref(space),
                make_shared_from_ref(space),
                make_shared_from_ref(space));

    arma::Mat<RT> matDefault = opDefault.weakForm()->asMatrix();
    arma::Mat<RT> matAccurate = opAccurate.weakForm()->asMatrix();

    const CT eps = std::numeric_limits<CT>::epsilon();
    BOOST_CHECK(check_arrays_are_close<RT>(
                    matDefault, matAccurate, 100 * eps));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(default_quadrature_order_is_sufficient_for_piecewise_quadratics,
                              BasisFunctionType, basis_function_types)
{
    typedef BasisFunctionType BFT;
    typedef BFT RT;
    typedef BFT KT;
    typedef typename Fiber::ScalarTraits<BFT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
                params, "meshes/sphere-ico-1.msh",
                false /* verbose */);

    PiecewisePolynomialContinuousScalarSpace<BFT> space(grid, 2);

    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    AccuracyOptions accuracyOptions;
    NumericalQuadratureStrategy<BFT, RT> defaultQuadStrategy(accuracyOptions);
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(5);
    accuracyOptions.doubleSingular.setRelativeQuadratureOrder(5);
    NumericalQuadratureStrategy<BFT, RT> accurateQuadStrategy(accuracyOptions);

    Context<BFT, RT> defaultContext(
                make_shared_from_ref(defaultQuadStrategy), assemblyOptions);
    Context<BFT, RT> accurateContext(
                make_shared_from_ref(accurateQuadStrategy), assemblyOptions);

    BoundaryOperator<BFT, RT> opDefault =
            constantBoundaryOperator<BFT, KT, RT>(
                make_shared_from_ref(defaultContext),
                make_shared_from_ref(space),
                make_shared_from_ref(space),
                make_shared_from_ref(space));
    BoundaryOperator<BFT, RT> opAccurate =
            constantBoundaryOperator<BFT, KT, RT>(
                make_shared_from_ref(accurateContext),
                make_shared_from_ref(space),
                make_shared_from_ref(space),
                make_shared_from_ref(space));

    arma::Mat<RT> matDefault = opDefault.weakForm()->asMatrix();
    arma::Mat<RT> matAccurate = opAccurate.weakForm()->asMatrix();

    const CT eps = std::numeric_limits<CT>::epsilon();
    BOOST_CHECK(check_arrays_are_close<RT>(
                    matDefault, matAccurate, 100 * eps));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(default_quadrature_order_is_sufficient_for_piecewise_cubics,
                              BasisFunctionType, basis_function_types)
{
    typedef BasisFunctionType BFT;
    typedef BFT RT;
    typedef BFT KT;
    typedef typename Fiber::ScalarTraits<BFT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
                params, "meshes/sphere-ico-1.msh",
                false /* verbose */);

    PiecewisePolynomialDiscontinuousScalarSpace<BFT> space(grid, 3);

    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    AccuracyOptions accuracyOptions;
    NumericalQuadratureStrategy<BFT, RT> defaultQuadStrategy(accuracyOptions);
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(5);
    accuracyOptions.doubleSingular.setRelativeQuadratureOrder(5);
    NumericalQuadratureStrategy<BFT, RT> accurateQuadStrategy(accuracyOptions);

    Context<BFT, RT> defaultContext(
                make_shared_from_ref(defaultQuadStrategy), assemblyOptions);
    Context<BFT, RT> accurateContext(
                make_shared_from_ref(accurateQuadStrategy), assemblyOptions);

    BoundaryOperator<BFT, RT> opDefault =
            constantBoundaryOperator<BFT, KT, RT>(
                make_shared_from_ref(defaultContext),
                make_shared_from_ref(space),
                make_shared_from_ref(space),
                make_shared_from_ref(space));
    BoundaryOperator<BFT, RT> opAccurate =
            constantBoundaryOperator<BFT, KT, RT>(
                make_shared_from_ref(accurateContext),
                make_shared_from_ref(space),
                make_shared_from_ref(space),
                make_shared_from_ref(space));

    arma::Mat<RT> matDefault = opDefault.weakForm()->asMatrix();
    arma::Mat<RT> matAccurate = opAccurate.weakForm()->asMatrix();

    const CT eps = std::numeric_limits<CT>::epsilon();
    BOOST_CHECK(check_arrays_are_close<RT>(
                    matDefault, matAccurate, 100 * eps));
}

BOOST_AUTO_TEST_SUITE_END()
