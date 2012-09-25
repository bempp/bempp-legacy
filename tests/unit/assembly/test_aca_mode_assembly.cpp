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
#include "assembly/discrete_boundary_operator.hpp"
#include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "grid/grid_factory.hpp"
#include "space/piecewise_constant_scalar_space.hpp"
#include "space/piecewise_linear_continuous_scalar_space.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/type_traits/is_complex.hpp>
#include "grid/grid.hpp"

using namespace Bempp;

// Tests

BOOST_AUTO_TEST_SUITE(AcaAssembly)

BOOST_AUTO_TEST_CASE_TEMPLATE(aca_of_assembled_operator_agrees_with_dense_assembly_for_614_element_mesh,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../examples/meshes/sphere-614.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseConstants(
        new PiecewiseConstantScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));
    pwiseConstants->assignDofs();
    pwiseLinears->assignDofs();

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(1);
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
            laplace3dDoubleLayerBoundaryOperator<BFT, RT>(
                contextAca, pwiseLinears, pwiseLinears, pwiseConstants);
    arma::Mat<RT> weakFormAca = opAca.weakForm()->asMatrix();

    BOOST_CHECK(check_arrays_are_close<ValueType>(
                    weakFormDense, weakFormAca, 2. * acaOptions.eps));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(aca_of_disassembled_operator_agrees_with_dense_assembly_for_614_element_mesh,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../examples/meshes/sphere-614.msh", false /* verbose */);

    shared_ptr<Space<BFT> > pwiseConstants(
        new PiecewiseConstantScalarSpace<BFT>(grid));
    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));
    pwiseConstants->assignDofs();
    pwiseLinears->assignDofs();

    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(1);
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
    acaOptions.globalAssemblyBeforeCompression = false;
    assemblyOptionsAca.switchToAcaMode(acaOptions);
    shared_ptr<Context<BFT, RT> > contextAca(
        new Context<BFT, RT>(quadStrategy, assemblyOptionsAca));

    BoundaryOperator<BFT, RT> opAca =
            laplace3dDoubleLayerBoundaryOperator<BFT, RT>(
                contextAca, pwiseLinears, pwiseLinears, pwiseConstants);
    arma::Mat<RT> weakFormAca = opAca.weakForm()->asMatrix();

    BOOST_CHECK(check_arrays_are_close<ValueType>(
                    weakFormDense, weakFormAca, 2. * acaOptions.eps));
}

BOOST_AUTO_TEST_SUITE_END()


