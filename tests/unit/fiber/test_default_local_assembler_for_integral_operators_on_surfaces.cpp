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

#include "../type_template.hpp"
#include "../check_arrays_are_close.hpp"

#include "assembly/context.hpp"
#include "assembly/general_elementary_singular_integral_operator_imp.hpp"
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "common/scalar_traits.hpp"
#include "fiber/geometrical_data.hpp"
#include "fiber/local_assembler_for_operators.hpp"
#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"
#include "grid/grid_view.hpp"
#include "space/piecewise_constant_scalar_space.hpp"
#include "space/piecewise_linear_continuous_scalar_space.hpp"

#include "grid/entity.hpp"
#include "grid/mapper.hpp"

#include <algorithm>
#include "common/armadillo_fwd.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/version.hpp>
#include <complex>

using namespace Bempp;

const int N_ELEMENTS_X = 2, N_ELEMENTS_Y = 3;

/** \brief Fixture class. */
template <typename BFT, typename RT>
class DefaultLocalAssemblerForIntegralOperatorsOnSurfacesManager
{
public:
    typedef typename ScalarTraits<RT>::RealType CT;
    typedef PiecewiseConstantScalarSpace<BFT> PiecewiseConstantSpace;
    typedef PiecewiseLinearContinuousScalarSpace<BFT> PiecewiseLinearSpace;
    typedef ElementaryIntegralOperator<BFT, CT, RT> Operator;
    typedef NumericalQuadratureStrategy<BFT, RT> QuadratureStrategy;
    typedef Fiber::RawGridGeometry<CT> RawGridGeometry;

    DefaultLocalAssemblerForIntegralOperatorsOnSurfacesManager(
            bool cacheSingularIntegrals)
    {
        // Create a Bempp grid
        shared_ptr<Grid> grid = createGrid();

        // Create context
        Fiber::AccuracyOptions options;
        options.doubleRegular.setRelativeQuadratureOrder(1);
        quadStrategy.reset(new QuadratureStrategy);

        AssemblyOptions assemblyOptions;
        assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
        assemblyOptions.enableSingularIntegralCaching(cacheSingularIntegrals);
        Context<BFT, RT> context(quadStrategy, assemblyOptions);

        // These important thing is that the domain and dualToRange spaces are
        // different
        piecewiseConstantSpace.reset(new PiecewiseConstantSpace(grid));
        piecewiseLinearSpace.reset(new PiecewiseLinearSpace(grid));

        BoundaryOperator<BFT, RT> bop =
                laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                    make_shared_from_ref(context),
                    piecewiseConstantSpace,
                    piecewiseLinearSpace,
                    piecewiseLinearSpace,
                    "SLP");
        op = boost::dynamic_pointer_cast<const Operator>(bop.abstractOperator());

        // Construct local assembler

        assembler = op->makeAssembler(*quadStrategy, assemblyOptions);
    }

private:
    shared_ptr<Grid> createGrid() {
        GridParameters params;
        params.topology = GridParameters::TRIANGULAR;

        const int dimGrid = 2;
        typedef double ctype;
        arma::Col<double> lowerLeft(dimGrid);
        arma::Col<double> upperRight(dimGrid);
        arma::Col<unsigned int> nElements(dimGrid);
        lowerLeft.fill(0);
        upperRight.fill(1);
        nElements(0) = N_ELEMENTS_X;
        nElements(1) = N_ELEMENTS_Y;

        return GridFactory::createStructuredGrid(
                    params, lowerLeft, upperRight, nElements);
    }

public:
    shared_ptr<PiecewiseConstantSpace> piecewiseConstantSpace;
    shared_ptr<PiecewiseLinearSpace> piecewiseLinearSpace;
    shared_ptr<const Operator> op;

    shared_ptr<QuadratureStrategy> quadStrategy;
    std::auto_ptr<typename Operator::LocalAssembler> assembler;
};

// Tests

BOOST_AUTO_TEST_SUITE(DefaultLocalAssemblerForIntegralOperatorsOnSurfaces)

template <typename ResultType>
void
both_variants_of_evaluateLocalWeakForms_agree_for_callVariant_TEST_TRIAL_and_cacheSingularIntegrals(
        bool cacheSingularIntegrals)
{
    DefaultLocalAssemblerForIntegralOperatorsOnSurfacesManager<
            typename ScalarTraits<ResultType>::RealType, ResultType> mgr(
                cacheSingularIntegrals);

    const int elementCount = N_ELEMENTS_X * N_ELEMENTS_Y * 2;
    const int testIndexCount = elementCount, trialIndexCount = elementCount;
    std::vector<int> testIndices(testIndexCount);
    for (int i = 0; i < testIndexCount; ++i)
        testIndices[i] = i;
    std::vector<int> trialIndices(trialIndexCount);
    for (int i = 0; i < trialIndexCount; ++i)
        trialIndices[i] = i;

    Fiber::_2dArray<arma::Mat<ResultType> > resultVariant2;
    mgr.assembler->evaluateLocalWeakForms(testIndices, trialIndices, resultVariant2);

    // Gather successive columns
    Fiber::_2dArray<arma::Mat<ResultType> > resultVariant1(
                testIndexCount, trialIndexCount);
    std::vector<arma::Mat<ResultType> > colResult;
    for (int trialI = 0; trialI < trialIndexCount; ++trialI)
    {
        mgr.assembler->evaluateLocalWeakForms(Fiber::TEST_TRIAL, testIndices,
                                              trialIndices[trialI],
                                              Fiber::ALL_DOFS, colResult);
        for (int testI = 0; testI < testIndexCount; ++testI)
            resultVariant1(testI, trialI) = colResult[testI];
    }

    BOOST_CHECK(check_arrays_are_close<ResultType>(
                    resultVariant1, resultVariant2, 1e-6));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(
        both_variants_of_evaluateLocalWeakForms_agree_for_callVariant_TEST_TRIAL_and_cacheSingularIntegrals_false,
        ResultType, result_types)
{
    both_variants_of_evaluateLocalWeakForms_agree_for_callVariant_TEST_TRIAL_and_cacheSingularIntegrals<
            ResultType>(false);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(
        both_variants_of_evaluateLocalWeakForms_agree_for_callVariant_TEST_TRIAL_and_cacheSingularIntegrals_true,
        ResultType, result_types)
{
    both_variants_of_evaluateLocalWeakForms_agree_for_callVariant_TEST_TRIAL_and_cacheSingularIntegrals<
            ResultType>(true);
}

template <typename ResultType>
void
both_variants_of_evaluateLocalWeakForms_agree_for_callVariant_TRIAL_TEST_and_cacheSingularIntegrals(
        bool cacheSingularIntegrals)
{
    DefaultLocalAssemblerForIntegralOperatorsOnSurfacesManager<
            typename ScalarTraits<ResultType>::RealType, ResultType> mgr(
                cacheSingularIntegrals);

    const int elementCount = N_ELEMENTS_X * N_ELEMENTS_Y * 2;
    const int testIndexCount = elementCount, trialIndexCount = elementCount;
    std::vector<int> testIndices(testIndexCount);
    for (int i = 0; i < testIndexCount; ++i)
        testIndices[i] = i;
    std::vector<int> trialIndices(trialIndexCount);
    for (int i = 0; i < trialIndexCount; ++i)
        trialIndices[i] = i;

    Fiber::_2dArray<arma::Mat<ResultType> > resultVariant2;
    mgr.assembler->evaluateLocalWeakForms(testIndices, trialIndices, resultVariant2);

    // Gather successive rows
    Fiber::_2dArray<arma::Mat<ResultType> > resultVariant1(
                testIndexCount, trialIndexCount);
    std::vector<arma::Mat<ResultType> > rowResult;
    for (int testI = 0; testI < testIndexCount; ++testI)
    {
        mgr.assembler->evaluateLocalWeakForms(Fiber::TRIAL_TEST, trialIndices,
                                              testIndices[testI],
                                              Fiber::ALL_DOFS, rowResult);
        for (int trialI = 0; trialI < trialIndexCount; ++trialI)
            resultVariant1(testI, trialI) = rowResult[trialI];
    }

    BOOST_CHECK(check_arrays_are_close<ResultType>(
                    resultVariant1, resultVariant2, 1e-6));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(
        both_variants_of_evaluateLocalWeakForms_agree_for_callVariant_TRIAL_TEST_and_cacheSingularIntegrals_false,
        ResultType, result_types)
{
    both_variants_of_evaluateLocalWeakForms_agree_for_callVariant_TRIAL_TEST_and_cacheSingularIntegrals<
            ResultType>(false);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(
        both_variants_of_evaluateLocalWeakForms_agree_for_callVariant_TRIAL_TEST_and_cacheSingularIntegrals_true,
        ResultType, result_types)
{
    both_variants_of_evaluateLocalWeakForms_agree_for_callVariant_TRIAL_TEST_and_cacheSingularIntegrals<
            ResultType>(true);
}

template <typename ResultType>
void
evaluateLocalWeakForms_for_ALL_DOFS_and_single_dof_agree_for_callVariant_TEST_TRIAL_and_cacheSingularIntegrals(
        bool cacheSingularIntegrals)
{
    DefaultLocalAssemblerForIntegralOperatorsOnSurfacesManager<
            typename ScalarTraits<ResultType>::RealType, ResultType> mgr(
                cacheSingularIntegrals);

    const int elementCount = N_ELEMENTS_X * N_ELEMENTS_Y * 2;
    std::vector<int> elementIndicesA(elementCount);
    for (int i = 0; i < elementCount; ++i)
        elementIndicesA[i] = i;

    int elementIndexB = 2;

    std::vector<arma::Mat<ResultType> > completeResult;
    mgr.assembler->evaluateLocalWeakForms(Fiber::TEST_TRIAL, elementIndicesA,
                                          elementIndexB,
                                          Fiber::ALL_DOFS, completeResult);
    int elementBDofCount = completeResult[0].n_cols;

    // Gather successive rows
    std::vector<arma::Mat<ResultType> > resultForSingleDof;
    std::vector<arma::Mat<ResultType> > expected(elementCount);
    for (int i = 0; i < elementCount; ++i)
        expected[i].set_size(completeResult[i].n_rows, completeResult[i].n_cols);

    for (int dof = 0; dof < elementBDofCount; ++dof)
    {
        mgr.assembler->evaluateLocalWeakForms(Fiber::TEST_TRIAL, elementIndicesA,
                                              elementIndexB,
                                              dof, resultForSingleDof);
        for (int i = 0; i < elementCount; ++i)
            expected[i].col(dof) = resultForSingleDof[i];
    }

    BOOST_CHECK(check_arrays_are_close<ResultType>(
                    completeResult, expected, 1e-6));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(
        evaluateLocalWeakForms_for_ALL_DOFS_and_single_dof_agree_for_callVariant_TEST_TRIAL_and_cacheSingularIntegrals_false,
        ResultType, result_types)
{
    evaluateLocalWeakForms_for_ALL_DOFS_and_single_dof_agree_for_callVariant_TEST_TRIAL_and_cacheSingularIntegrals<
            ResultType>(false);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(
        evaluateLocalWeakForms_for_ALL_DOFS_and_single_dof_agree_for_callVariant_TEST_TRIAL_and_cacheSingularIntegrals_true,
        ResultType, result_types)
{
    evaluateLocalWeakForms_for_ALL_DOFS_and_single_dof_agree_for_callVariant_TEST_TRIAL_and_cacheSingularIntegrals<
            ResultType>(true);
}

template <typename ResultType>
void
evaluateLocalWeakForms_for_ALL_DOFS_and_single_dof_agree_for_callVariant_TRIAL_TEST_and_cacheSingularIntegrals(
        bool cacheSingularIntegrals)
{
    DefaultLocalAssemblerForIntegralOperatorsOnSurfacesManager<
            typename ScalarTraits<ResultType>::RealType, ResultType> mgr(
                cacheSingularIntegrals);

    const int elementCount = N_ELEMENTS_X * N_ELEMENTS_Y * 2;
    std::vector<int> elementIndicesA(elementCount);
    for (int i = 0; i < elementCount; ++i)
        elementIndicesA[i] = i;

    int elementIndexB = 2;

    std::vector<arma::Mat<ResultType> > completeResult;
    mgr.assembler->evaluateLocalWeakForms(Fiber::TRIAL_TEST, elementIndicesA,
                                          elementIndexB,
                                          Fiber::ALL_DOFS, completeResult);
    int elementBDofCount = completeResult[0].n_rows;

    // Gather successive rows
    std::vector<arma::Mat<ResultType> > resultForSingleDof;
    std::vector<arma::Mat<ResultType> > expected(elementCount);
    for (int i = 0; i < elementCount; ++i)
        expected[i].set_size(completeResult[i].n_rows, completeResult[i].n_cols);

    for (int dof = 0; dof < elementBDofCount; ++dof)
    {
        mgr.assembler->evaluateLocalWeakForms(Fiber::TRIAL_TEST, elementIndicesA,
                                              elementIndexB,
                                              dof, resultForSingleDof);
        for (int i = 0; i < elementCount; ++i)
            expected[i].row(dof) = resultForSingleDof[i];
    }

    BOOST_CHECK(check_arrays_are_close<ResultType>(
                    completeResult, expected, 1e-6));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(
        evaluateLocalWeakForms_for_ALL_DOFS_and_single_dof_agree_for_callVariant_TRIAL_TEST_and_cacheSingularIntegrals_false,
        ResultType, result_types)
{
    evaluateLocalWeakForms_for_ALL_DOFS_and_single_dof_agree_for_callVariant_TRIAL_TEST_and_cacheSingularIntegrals<
            ResultType>(false);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(
        evaluateLocalWeakForms_for_ALL_DOFS_and_single_dof_agree_for_callVariant_TRIAL_TEST_and_cacheSingularIntegrals_true,
        ResultType, result_types)
{
    evaluateLocalWeakForms_for_ALL_DOFS_and_single_dof_agree_for_callVariant_TRIAL_TEST_and_cacheSingularIntegrals<
            ResultType>(true);
}

template <typename ResultType>
void
evaluateLocalWeakForms_with_and_without_singular_integral_caching_gives_same_results(
        bool cacheSingularIntegrals)
{
    const int elementCount = N_ELEMENTS_X * N_ELEMENTS_Y * 2;
    const int testIndexCount = elementCount, trialIndexCount = elementCount;
    std::vector<int> testIndices(testIndexCount);
    for (int i = 0; i < testIndexCount; ++i)
        testIndices[i] = i;
    std::vector<int> trialIndices(trialIndexCount);
    for (int i = 0; i < trialIndexCount; ++i)
        trialIndices[i] = i;

    Fiber::_2dArray<arma::Mat<ResultType> > resultWithCaching;
    Fiber::_2dArray<arma::Mat<ResultType> > resultWithoutCaching;

    {
        DefaultLocalAssemblerForIntegralOperatorsOnSurfacesManager<
                typename ScalarTraits<ResultType>::RealType, ResultType> mgr(
                    true);
        mgr.assembler->evaluateLocalWeakForms(testIndices, trialIndices,
                                              resultWithCaching);
    }
    {
        DefaultLocalAssemblerForIntegralOperatorsOnSurfacesManager<
            typename ScalarTraits<ResultType>::RealType, ResultType> mgr(
                false);
        mgr.assembler->evaluateLocalWeakForms(testIndices, trialIndices,
                                              resultWithoutCaching);
    }

        BOOST_CHECK(check_arrays_are_close<ResultType>(
                    resultWithCaching, resultWithoutCaching, 1e-6));
}

BOOST_AUTO_TEST_SUITE_END()
