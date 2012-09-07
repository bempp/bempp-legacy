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

#ifndef bempp_laplace_3d_dirichlet_fixture
#define bempp_laplace_3d_dirichlet_fixture

#include "common/common.hpp"

#include "assembly/boundary_operator.hpp"
#include "assembly/grid_function.hpp"

#include "common/armadillo_fwd.hpp"
#include "common/scalar_traits.hpp"

#include "grid/grid.hpp"

#include <memory>

namespace Bempp
{

template <typename ValueType_>
struct UnitFunctor
{
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    int argumentDimension() const { return 3; }
    int resultDimension() const { return 1; }

    void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        result(0) = 1.;
    }
};

enum SpaceType { 
    PIECEWISE_CONSTANTS,
    PIECEWISE_LINEARS
};

template <typename BFT, typename RT>
struct Laplace3dDirichletFixture
{
    Laplace3dDirichletFixture(
        SpaceType dirichletDataDomain = PIECEWISE_LINEARS,
        SpaceType neumannDataDomain = PIECEWISE_CONSTANTS,
        SpaceType range = PIECEWISE_LINEARS,
        SpaceType dualToRange = PIECEWISE_CONSTANTS);

    shared_ptr<Grid> grid;
    BoundaryOperator<BFT, RT> lhsOp;
    GridFunction<BFT, RT> rhs;
};

} // namespace Bempp

#endif
