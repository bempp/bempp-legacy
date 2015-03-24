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
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"

#include "common/boost_make_shared_fwd.hpp"

#include "grid/grid_factory.hpp"
#include "grid/grid.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"

#include <algorithm>
#include "common/eigen_support.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/version.hpp>

// Tests

using namespace Bempp;

BOOST_AUTO_TEST_SUITE(Laplace3dSingleLayerBoundaryOperator)

BOOST_AUTO_TEST_CASE_TEMPLATE(symmetric_matches_nonsymmetric_in_aca_mode,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    if (boost::is_same<RT, std::complex<float> >::value) {
        // The AHMED support for single-precision complex symmetric matrices
        // is broken
        BOOST_CHECK(true);
        return;
    }

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
                params, "../../meshes/sphere-h-0.4.msh",
                false /* verbose */);

    PiecewiseLinearContinuousScalarSpace<BFT> pwiseLinears(grid);
    PiecewiseConstantScalarSpace<BFT> pwiseConstants(grid);

    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    AcaOptions acaOptions;
    acaOptions.minimumBlockSize = 4;
    assemblyOptions.switchToAcaMode(acaOptions);
    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(4);
    accuracyOptions.doubleSingular.setRelativeQuadratureOrder(2);
    NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);

    Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);
 
    BoundaryOperator<BFT, RT> opNonsymmetric =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseConstants),
                make_shared_from_ref(pwiseLinears),
                "", NO_SYMMETRY);
    BoundaryOperator<BFT, RT> opSymmetric =
            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseConstants),
                make_shared_from_ref(pwiseLinears),
                "", SYMMETRIC);
 
#   ifdef WITH_AHMED
       arma::Mat<RT> matNonsymmetric = opNonsymmetric.weakForm()->asMatrix();
       arma::Mat<RT> matSymmetric = opSymmetric.weakForm()->asMatrix();
    
       BOOST_CHECK(check_arrays_are_close<RT>(
                       matNonsymmetric, matSymmetric, 2 * acaOptions.eps));
#   else
       BOOST_CHECK(true);
#   endif
}

BOOST_AUTO_TEST_SUITE_END()
