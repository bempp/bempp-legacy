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

#include "../check_arrays_are_close.hpp"
#include "../type_template.hpp"

#include <algorithm>
#include "common/armadillo_fwd.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/version.hpp>

// Tests

using namespace Bempp;

BOOST_AUTO_TEST_SUITE(SurfaceCurl3dFunctor)

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_values_works,
                              ValueType, basis_function_types)
{
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
    typedef Fiber::SurfaceCurl3dElementaryFunctor<CoordinateType> ElementaryFunctor;
    typedef Fiber::ElementaryBasisTransformationFunctorWrapper<ElementaryFunctor> Functor;
    Functor functor;
    size_t basisDeps = 0, geomDeps = 0;
    functor.addDependencies(basisDeps, geomDeps);

    arma::Mat<CoordinateType> points(2, 4);
    srand(1);
    points.randu();

    typedef Fiber::PiecewiseLinearContinuousScalarBasis<3, ValueType> Basis;
    Basis basis;
    Fiber::BasisData<ValueType> basisData;
    basis.evaluate(basisDeps, points, Fiber::ALL_DOFS, basisData);

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
                params, "meshes/simple_mesh_2_elements.msh",
                false /* verbose */);
    std::unique_ptr<GridView> view = grid->leafView();
    std::unique_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    it->next();
    const Entity<0>& element = it->entity();
    const Geometry& geo = element.geometry();

    Fiber::GeometricalData<CoordinateType> geomData;
    geo.getData(geomDeps, points, geomData);

    Fiber::DefaultCollectionOfBasisTransformations<Functor> transformations(functor);

    Fiber::CollectionOf3dArrays<ValueType> result;
    transformations.evaluate(basisData, geomData, result);

    Fiber::_3dArray<ValueType> expected(functor.resultDimension(0),
                                        basisData.functionCount(),
                                        basisData.pointCount());
    for (size_t p = 0; p < basisData.pointCount(); ++p) {
        // component, function, point
        expected(0, 0, p) = 1.;
        expected(1, 0, p) = 0.;
        expected(2, 0, p) = 0.;
        expected(0, 1, p) = -5./7.;
        expected(1, 1, p) = 10./7.;
        expected(2, 1, p) = 0.;
        expected(0, 2, p) = -2./7.;
        expected(1, 2, p) = -10./7.;
        expected(2, 2, p) = 0.;
    }

    BOOST_CHECK(check_arrays_are_close<ValueType>(
                    expected, result[0], 1e-13));
}

BOOST_AUTO_TEST_SUITE_END()
