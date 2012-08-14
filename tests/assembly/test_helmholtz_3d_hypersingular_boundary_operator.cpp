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
#include "assembly/helmholtz_3d_single_layer_boundary_operator.hpp"
#include "assembly/helmholtz_3d_hypersingular_boundary_operator.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"

#include "common/boost_make_shared_fwd.hpp"

#include "fiber/basis_transformation_functor_wrappers.hpp"
#include "fiber/default_collection_of_basis_transformations.hpp"
#include "fiber/geometrical_data.hpp"
#include "fiber/piecewise_linear_continuous_scalar_basis.hpp"
#include "fiber/surface_curl_3d_functor.hpp"

#include "grid/grid_factory.hpp"
#include "grid/grid.hpp"
#include "grid/grid_view.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/entity.hpp"
#include "grid/geometry.hpp"

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

BOOST_AUTO_TEST_SUITE(Helmholtz3dHypersingularBoundaryOperator)

BOOST_AUTO_TEST_CASE_TEMPLATE(works, BasisFunctionType, basis_function_types)
{
    typedef BasisFunctionType BFT;
    typedef typename Fiber::ScalarTraits<BFT>::ComplexType RT;
    typedef typename Fiber::ScalarTraits<BFT>::RealType CT;
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    std::auto_ptr<Grid> grid = GridFactory::importGmshGrid(
                params, std::string("../examples/meshes/two_disjoint_triangles.msh"));

    PiecewiseLinearContinuousScalarSpace<BFT> pwiseLinears(*grid);
    PiecewiseConstantScalarSpace<BFT> pwiseConstants(*grid);
    pwiseLinears.assignDofs();
    pwiseConstants.assignDofs();

    AssemblyOptions assemblyOptions;
    AccuracyOptions accuracyOptions;
    accuracyOptions.doubleRegular.orderIncrement = 5;
    NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);

    Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

    const RT waveNumber(1.23, 0.31);

    BoundaryOperator<BFT, RT> slpOpConstants = Bempp::helmholtz3dSingleLayerBoundaryOperator<BFT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseConstants),
                make_shared_from_ref(pwiseConstants),
                make_shared_from_ref(pwiseConstants),
                waveNumber);
    BoundaryOperator<BFT, RT> slpOpLinears = helmholtz3dSingleLayerBoundaryOperator<BFT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                waveNumber);
    BoundaryOperator<BFT, RT> hypOp = helmholtz3dHypersingularBoundaryOperator<BFT>(
                make_shared_from_ref(context),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                make_shared_from_ref(pwiseLinears),
                waveNumber);

    // Get the matrix repr. of the hypersingular operator
    arma::Mat<RT> hypMat = hypOp.weakForm()->asMatrix();

    // Construct the expected hypersingular operator matrix. For this, we need:

    // * the surface curls of all basis functions (which are constant)
    typedef Fiber::SurfaceCurl3dElementaryFunctor<CT> ElementaryFunctor;
    typedef Fiber::ElementaryBasisTransformationFunctorWrapper<ElementaryFunctor> Functor;
    Functor functor;
    size_t basisDeps = 0, geomDeps = 0;
    functor.addDependencies(basisDeps, geomDeps);

    arma::Mat<CT> points(2, 1);
    points.fill(0.);

    typedef Fiber::PiecewiseLinearContinuousScalarBasis<3, BFT> Basis;
    Basis basis;
    Fiber::BasisData<BFT> basisData;
    basis.evaluate(basisDeps, points, Fiber::ALL_DOFS, basisData);

    Fiber::GeometricalData<CT> geomData[2];
    std::auto_ptr<GridView> view = grid->leafView();
    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    it->entity().geometry().getData(geomDeps, points, geomData[0]);
    it->next();
    it->entity().geometry().getData(geomDeps, points, geomData[1]);

    Fiber::DefaultCollectionOfBasisTransformations<Functor> transformations(functor);

    Fiber::CollectionOf3dArrays<BFT> surfaceCurl[2];
    transformations.evaluate(basisData, geomData[0], surfaceCurl[0]);
    transformations.evaluate(basisData, geomData[1], surfaceCurl[1]);

    // * the single-layer potential matrix for constant basis functions
    arma::Mat<RT> slpMatConstants = slpOpConstants.weakForm()->asMatrix();
    // * the single-layer potential matrix for linear basis functions
    arma::Mat<RT> slpMatLinears = slpOpLinears.weakForm()->asMatrix();

    arma::Mat<RT> expectedHypMat(6, 6);
    for (size_t testElement = 0; testElement < 2; ++testElement)
        for (size_t trialElement = 0; trialElement < 2; ++trialElement)
            for (size_t r = 0; r < 3; ++r)
                for (size_t c = 0; c < 3; ++c) {
                    RT curlMultiplier = 0.;
                    for (size_t dim = 0; dim < 3; ++dim)
                        curlMultiplier += surfaceCurl[testElement][0](dim, r, 0) *
                                surfaceCurl[trialElement][0](dim, c, 0);
                    RT valueMultiplier = 0.;
                    for (size_t dim = 0; dim < 3; ++dim)
                        valueMultiplier += geomData[testElement].normals(dim, 0) * 
                            geomData[trialElement].normals(dim, 0);
                    expectedHypMat(3 * testElement + r, 3 * trialElement + c) =
                            curlMultiplier * slpMatConstants(testElement, trialElement) -
                            waveNumber * waveNumber * valueMultiplier *
                            slpMatLinears(3 * testElement + r, 3 * trialElement + c);
                }

    BOOST_CHECK(check_arrays_are_close<RT>(expectedHypMat, hypMat, 1e-6));
}

BOOST_AUTO_TEST_SUITE_END()
