// Copyright (C) 2011-2012 by the BEM++ Authors
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

#include "piecewise_linear_scalar_space.hpp"

#include "space_helper.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../grid/entity.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid.hpp"

#include <stdexcept>
#include <iostream>

namespace Bempp
{

template <typename BasisFunctionType>
PiecewiseLinearScalarSpace<BasisFunctionType>::
PiecewiseLinearScalarSpace(const shared_ptr<const Grid>& grid) :
    ScalarSpace<BasisFunctionType>(grid)
{
}

template <typename BasisFunctionType>
PiecewiseLinearScalarSpace<BasisFunctionType>::
~PiecewiseLinearScalarSpace()
{
}

template <typename BasisFunctionType>
int PiecewiseLinearScalarSpace<BasisFunctionType>::domainDimension() const
{
    return this->grid()->dim();
}

template <typename BasisFunctionType>
int PiecewiseLinearScalarSpace<BasisFunctionType>::codomainDimension() const
{
    return 1;
}

template <typename BasisFunctionType>
const Fiber::Shapeset<BasisFunctionType>&
PiecewiseLinearScalarSpace<BasisFunctionType>::shapeset(
        const Entity<0>& element) const
{
    switch (elementVariant(element))
    {
    case 3:
        return m_triangleShapeset;
    case 4:
        return m_quadrilateralShapeset;
    case 2:
        return m_lineShapeset;
    default:
        throw std::logic_error("PiecewiseLinearScalarSpace::shapeset(): "
                               "invalid element variant, this shouldn't happen!");
    }
}

template <typename BasisFunctionType>
ElementVariant PiecewiseLinearScalarSpace<BasisFunctionType>::elementVariant(
        const Entity<0>& element) const
{
    GeometryType type = element.type();
    if (type.isLine())
        return 2;
    else if (type.isTriangle())
        return 3;
    else if (type.isQuadrilateral())
        return 4;
    else
        throw std::runtime_error("PiecewiseLinearScalarSpace::"
                                 "elementVariant(): invalid geometry type, "
                                 "this shouldn't happen!");
}

template <typename BasisFunctionType>
void PiecewiseLinearScalarSpace<BasisFunctionType>::setElementVariant(
        const Entity<0>& element, ElementVariant variant)
{
    if (variant != elementVariant(element))
        // for this space, the element variants are unmodifiable,
        throw std::runtime_error("PiecewiseLinearScalarSpace::"
                                 "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
void PiecewiseLinearScalarSpace<BasisFunctionType>::
getGlobalDofInterpolationPoints(arma::Mat<CoordinateType>& points) const
{
    SpaceHelper<BasisFunctionType>::
            getGlobalDofInterpolationPoints_defaultImplementation(
                *this, points);
}

template <typename BasisFunctionType>
void PiecewiseLinearScalarSpace<BasisFunctionType>::
getNormalsAtGlobalDofInterpolationPoints(arma::Mat<CoordinateType>& normals) const
{
    SpaceHelper<BasisFunctionType>::
            getNormalsAtGlobalDofInterpolationPoints_defaultImplementation(
                *this, normals);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(PiecewiseLinearScalarSpace);

} // namespace Bempp
