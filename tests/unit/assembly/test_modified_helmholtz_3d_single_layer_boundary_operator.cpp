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
#include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
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
#include <complex>

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

using namespace Bempp;

BOOST_AUTO_TEST_SUITE(ModifiedHelmholtz3dSingleLayerBoundaryOperator)

BOOST_AUTO_TEST_CASE_TEMPLATE(interpolated_matches_noniterpolated,
                              BasisFunctionType, basis_function_types)
{
    typedef BasisFunctionType BFT;
    typedef typename Fiber::ScalarTraits<BFT>::ComplexType RT;
    typedef typename Fiber::ScalarTraits<BFT>::ComplexType KT;
    typedef typename Fiber::ScalarTraits<BFT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
                params, "meshes/two_disjoint_triangles.msh",
                false /* verbose */);

    PiecewiseLinearContinuousScalarSpace<BFT> pwiseLinears(grid);
    PiecewiseConstantScalarSpace<BFT> pwiseConstants(grid);

    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.setAbsoluteQuadratureOrder(5);
    accuracyOptions.doubleSingular.setAbsoluteQuadratureOrder(5);
    NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);

    Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

    const KT waveNumber(3.23, 0.31);

    BoundaryOperator<BFT, RT> opNoninterpolated =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, KT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                waveNumber,
                "", NO_SYMMETRY,
                false);
    BoundaryOperator<BFT, RT> opInterpolated =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, KT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                waveNumber,
                "", NO_SYMMETRY,
                true);

    Matrix<RT> matNoninterpolated = opNoninterpolated.weakForm()->asMatrix();
    Matrix<RT> matInterpolated = opInterpolated.weakForm()->asMatrix();

    const CT eps = std::numeric_limits<CT>::epsilon();
    BOOST_CHECK(check_arrays_are_close<RT>(
                    matNoninterpolated, matInterpolated, 100 * eps));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(symmetric_matches_nonsymmetric_in_aca_mode,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename Fiber::ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

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

    const RT waveNumber = initWaveNumber<RT>();

    BoundaryOperator<BFT, RT> opNonsymmetric =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, RT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseConstants),
                make_shared_from_ref(pwiseLinears),
                waveNumber,
                "", NO_SYMMETRY);
    BoundaryOperator<BFT, RT> opSymmetric =
            modifiedHelmholtz3dSingleLayerBoundaryOperator<BFT, RT, RT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseConstants),
                make_shared_from_ref(pwiseLinears),
                waveNumber,
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
