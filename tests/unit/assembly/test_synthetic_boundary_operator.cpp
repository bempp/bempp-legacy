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
#include "../check_arrays_are_close.hpp"

#include "assembly/context.hpp"
#include "assembly/discrete_boundary_operator.hpp"
#include "assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_hypersingular_boundary_operator.hpp"
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/helmholtz_3d_hypersingular_boundary_operator.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "grid/grid_factory.hpp"
#include "space/piecewise_constant_scalar_space.hpp"
#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_linear_discontinuous_scalar_space.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/type_traits/is_complex.hpp>
#include "grid/grid.hpp"

using namespace Bempp;

// Tests

BOOST_AUTO_TEST_SUITE(SyntheticBoundaryOperator)

BOOST_AUTO_TEST_CASE_TEMPLATE(aca_of_synthetic_single_layer_operator_agrees_with_dense_assembly_for_pwise_constants,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../examples/meshes/sphere-h-0.4.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseConstants(
        new PiecewiseConstantScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseDLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsDense;
    assemblyOptionsDense.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextDense(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsDense));

    BoundaryOperator<BFT, RT> opDense =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                contextDense, pwiseConstants, pwiseLinears, pwiseConstants);
    arma::Mat<RT> weakFormDense = opDense.weakForm()->asMatrix();

    AssemblyOptions assemblyOptionsAca;
    assemblyOptionsAca.setVerbosityLevel(VerbosityLevel::LOW);
    AcaOptions acaOptions;
    assemblyOptionsAca.switchToAcaMode(acaOptions);
    shared_ptr<Context<BFT, RT> > contextAca(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsAca));

    BoundaryOperator<BFT, RT> opAca =
            laplace3dSyntheticSingleLayerBoundaryOperator<BFT, RT>(
                contextAca, pwiseConstants, pwiseLinears, pwiseConstants,
                pwiseConstants, pwiseConstants);
    arma::Mat<RT> weakFormAca = opAca.weakForm()->asMatrix();

    BOOST_CHECK(check_arrays_are_close<ValueType>(
                    weakFormDense, weakFormAca, 2. * acaOptions.eps));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(aca_of_synthetic_single_layer_operator_agrees_with_dense_assembly_for_pwise_linears,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../examples/meshes/sphere-h-0.4.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseConstants(
        new PiecewiseConstantScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseDLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsDense;
    assemblyOptionsDense.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextDense(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsDense));

    BoundaryOperator<BFT, RT> opDense =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                contextDense, pwiseLinears, pwiseLinears, pwiseLinears);
    arma::Mat<RT> weakFormDense = opDense.weakForm()->asMatrix();

    AssemblyOptions assemblyOptionsAca;
    assemblyOptionsAca.setVerbosityLevel(VerbosityLevel::LOW);
    AcaOptions acaOptions;
    assemblyOptionsAca.switchToAcaMode(acaOptions);
    shared_ptr<Context<BFT, RT> > contextAca(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsAca));

    BoundaryOperator<BFT, RT> opAca =
            laplace3dSyntheticSingleLayerBoundaryOperator<BFT, RT>(
                contextAca, pwiseLinears, pwiseLinears, pwiseLinears,
                pwiseDLinears, pwiseDLinears);
    arma::Mat<RT> weakFormAca = opAca.weakForm()->asMatrix();

    BOOST_CHECK(check_arrays_are_close<ValueType>(
                    weakFormDense, weakFormAca, 2. * acaOptions.eps));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(aca_of_synthetic_double_layer_operator_agrees_with_dense_assembly,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../examples/meshes/sphere-h-0.4.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseConstants(
        new PiecewiseConstantScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseDLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsDense;
    assemblyOptionsDense.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextDense(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsDense));

    BoundaryOperator<BFT, RT> opDense =
            laplace3dDoubleLayerBoundaryOperator<BFT, RT>(
                contextDense, pwiseLinears, pwiseLinears, pwiseConstants);
    arma::Mat<RT> weakFormDense = opDense.weakForm()->asMatrix();

    AssemblyOptions assemblyOptionsAca;
    assemblyOptionsAca.setVerbosityLevel(VerbosityLevel::LOW);
    AcaOptions acaOptions;
    assemblyOptionsAca.switchToAcaMode(acaOptions);
    shared_ptr<Context<BFT, RT> > contextAca(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsAca));

    BoundaryOperator<BFT, RT> opAca =
            laplace3dSyntheticDoubleLayerBoundaryOperator<BFT, RT>(
                contextAca, pwiseLinears, pwiseLinears, pwiseConstants,
                pwiseDLinears, pwiseConstants);
    arma::Mat<RT> weakFormAca = opAca.weakForm()->asMatrix();

    BOOST_CHECK(check_arrays_are_close<ValueType>(
                    weakFormDense, weakFormAca, 2. * acaOptions.eps));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(aca_of_synthetic_adjoint_double_layer_operator_agrees_with_dense_assembly,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../examples/meshes/sphere-h-0.4.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseConstants(
        new PiecewiseConstantScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseDLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsDense;
    assemblyOptionsDense.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextDense(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsDense));

    BoundaryOperator<BFT, RT> opDense =
            laplace3dAdjointDoubleLayerBoundaryOperator<BFT, RT>(
                contextDense, pwiseConstants, pwiseConstants, pwiseLinears);
    arma::Mat<RT> weakFormDense = opDense.weakForm()->asMatrix();

    AssemblyOptions assemblyOptionsAca;
    assemblyOptionsAca.setVerbosityLevel(VerbosityLevel::LOW);
    AcaOptions acaOptions;
    assemblyOptionsAca.switchToAcaMode(acaOptions);
    shared_ptr<Context<BFT, RT> > contextAca(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsAca));

    BoundaryOperator<BFT, RT> opAca =
            laplace3dSyntheticAdjointDoubleLayerBoundaryOperator<BFT, RT>(
                contextAca, pwiseConstants, pwiseConstants, pwiseLinears,
                pwiseConstants, pwiseDLinears);
    arma::Mat<RT> weakFormAca = opAca.weakForm()->asMatrix();

    BOOST_CHECK(check_arrays_are_close<ValueType>(
                    weakFormDense, weakFormAca, 2. * acaOptions.eps));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(aca_of_synthetic_hypersingular_operator_agrees_with_dense_assembly_in_symmetric_case,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../examples/meshes/sphere-h-0.4.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseConstants(
        new PiecewiseConstantScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseDLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsDense;
    assemblyOptionsDense.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextDense(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsDense));

    BoundaryOperator<BFT, RT> opDense =
            laplace3dHypersingularBoundaryOperator<BFT, RT>(
                contextDense, pwiseLinears, pwiseConstants, pwiseLinears);
    arma::Mat<RT> weakFormDense = opDense.weakForm()->asMatrix();

    AssemblyOptions assemblyOptionsAca;
    assemblyOptionsAca.setVerbosityLevel(VerbosityLevel::LOW);
    AcaOptions acaOptions;
    assemblyOptionsAca.switchToAcaMode(acaOptions);
    shared_ptr<Context<BFT, RT> > contextAca(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsAca));

    BoundaryOperator<BFT, RT> opAca =
            laplace3dSyntheticHypersingularBoundaryOperator<BFT, RT>(
                contextAca, pwiseLinears, pwiseConstants, pwiseLinears,
                pwiseDLinears, pwiseDLinears);
    arma::Mat<RT> weakFormAca = opAca.weakForm()->asMatrix();

    BOOST_CHECK(check_arrays_are_close<ValueType>(
                    weakFormDense, weakFormAca, 2. * acaOptions.eps));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(aca_of_synthetic_hypersingular_operator_agrees_with_dense_assembly_in_asymmetric_case,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../examples/meshes/sphere-h-0.4.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseConstants(
        new PiecewiseConstantScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseDLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseDLinears2(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsDense;
    assemblyOptionsDense.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > contextDense(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsDense));

    BoundaryOperator<BFT, RT> opDense =
            laplace3dHypersingularBoundaryOperator<BFT, RT>(
                contextDense, pwiseLinears, pwiseConstants, pwiseLinears);
    arma::Mat<RT> weakFormDense = opDense.weakForm()->asMatrix();

    AssemblyOptions assemblyOptionsAca;
    assemblyOptionsAca.setVerbosityLevel(VerbosityLevel::LOW);
    AcaOptions acaOptions;
    assemblyOptionsAca.switchToAcaMode(acaOptions);
    shared_ptr<Context<BFT, RT> > contextAca(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsAca));

    // Internal domain different from dualToRange
    BoundaryOperator<BFT, RT> opAca =
            laplace3dSyntheticHypersingularBoundaryOperator<BFT, RT>(
                contextAca, pwiseLinears, pwiseConstants, pwiseLinears,
                pwiseDLinears, pwiseDLinears2);
    arma::Mat<RT> weakFormAca = opAca.weakForm()->asMatrix();

    BOOST_CHECK(check_arrays_are_close<ValueType>(
                    weakFormDense, weakFormAca, 2. * acaOptions.eps));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(aca_of_synthetic_hypersingular_operator_throws_for_continuous_internal_domain,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../examples/meshes/sphere-h-0.4.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseConstants(
        new PiecewiseConstantScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseDLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsAca;
    assemblyOptionsAca.setVerbosityLevel(VerbosityLevel::LOW);
    AcaOptions acaOptions;
    assemblyOptionsAca.switchToAcaMode(acaOptions);
    shared_ptr<Context<BFT, RT> > contextAca(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsAca));

    BoundaryOperator<BFT, RT> opAca =
            laplace3dSyntheticHypersingularBoundaryOperator<BFT, RT>(
                contextAca, pwiseLinears, pwiseConstants, pwiseLinears,
                pwiseLinears, pwiseDLinears);
    BOOST_CHECK_THROW(opAca.weakForm(), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(aca_of_synthetic_hypersingular_operator_throws_for_continuous_internal_dualToRange,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../examples/meshes/sphere-h-0.4.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseConstants(
        new PiecewiseConstantScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseDLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsAca;
    assemblyOptionsAca.setVerbosityLevel(VerbosityLevel::LOW);
    AcaOptions acaOptions;
    assemblyOptionsAca.switchToAcaMode(acaOptions);
    shared_ptr<Context<BFT, RT> > contextAca(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsAca));

    BoundaryOperator<BFT, RT> opAca =
            laplace3dSyntheticHypersingularBoundaryOperator<BFT, RT>(
                contextAca, pwiseLinears, pwiseConstants, pwiseLinears,
                pwiseDLinears, pwiseLinears);
    BOOST_CHECK_THROW(opAca.weakForm(), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(aca_of_synthetic_hypersingular_operator_throws_for_continuous_internal_domain_and_dualToRange,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../examples/meshes/sphere-h-0.4.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseConstants(
        new PiecewiseConstantScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseDLinears(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2);
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));

    AssemblyOptions assemblyOptionsAca;
    assemblyOptionsAca.setVerbosityLevel(VerbosityLevel::LOW);
    AcaOptions acaOptions;
    assemblyOptionsAca.switchToAcaMode(acaOptions);
    shared_ptr<Context<BFT, RT> > contextAca(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsAca));

    BoundaryOperator<BFT, RT> opAca =
            laplace3dSyntheticHypersingularBoundaryOperator<BFT, RT>(
                contextAca, pwiseLinears, pwiseConstants, pwiseLinears,
                pwiseLinears, pwiseLinears);
    BOOST_CHECK_THROW(opAca.weakForm(), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

#endif // WITH_AHMED
