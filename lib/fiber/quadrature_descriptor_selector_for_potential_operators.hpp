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

#ifndef fiber_quadrature_descriptor_selector_for_potential_operators_hpp
#define fiber_quadrature_descriptor_selector_for_potential_operators_hpp

#include "../common/common.hpp"

#include "../common/armadillo_fwd.hpp"
#include "scalar_traits.hpp"
#include "single_quadrature_descriptor.hpp"

namespace Fiber
{

template <typename BasisFunctionType> class Shapeset;

/** \ingroup quadrature
 *  \brief Quadrature descriptor selector used during the
 *  evaluation on potentials. */
template <typename BasisFunctionType>
class QuadratureDescriptorSelectorForPotentialOperators
{
public:
    /** \brief Type used to represent coordinates. */
    typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;
    /** \brief Destructor. */
    virtual ~QuadratureDescriptorSelectorForPotentialOperators() {}

    /** \brief Return the descriptor of the quadrature rule to be used
     *  in the evaluation of the contribution of element
     *  \p trialElementIndex at point \p point.
     *
     *  \param[in] point
     *    Coordinates of the point where the potential is evaluated.
     *  \param[in] trialElementIndex
     *    Index of the element whose contribution to the potential is evaluated.
     *  \param[in] nominalDistance
     *    This parameter is only relevant for quadrature selectors
     *    that vary the degree of accuracy of regular quadrature rules
     *    with point-element distance. If \p nominalDistance is
     *    non-negative, the quadrature rule should be chosen as if the
     *    distance between the evaluation point and the center of the
     *    element were equal to \p nominalDistance. Otherwise the
     *    selector should try to estimate the distance on its own and
     *    adjust the quadrature order accordingly. */
    virtual SingleQuadratureDescriptor quadratureDescriptor(
        const arma::Col<CoordinateType>& point, int trialElementIndex,
        CoordinateType nominalDistance) const = 0;

    /** \brief Return the descriptor of the quadrature rule used in
     *  the evaluation of a potential "far away" from the surface.
     *
     *  This function should return the quadrature description to be
     *  used in the evaluation of the contribution of an element with
     *  \p trialElementCornerCount vertices, endowed with the \p
     *  trialShapeset set of shape functions, at points lying "far
     *  away" from the element. */
    virtual SingleQuadratureDescriptor farFieldQuadratureDescriptor(
        const Shapeset<BasisFunctionType>& trialShapeset,
        int trialElementCornerCount) const = 0;
};

} // namespace Fiber

#endif
