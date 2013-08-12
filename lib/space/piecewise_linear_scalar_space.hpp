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

#ifndef bempp_piecewise_linear_scalar_space_hpp
#define bempp_piecewise_linear_scalar_space_hpp

#include "../common/common.hpp"

#include "../grid/grid_view.hpp"
#include "scalar_space.hpp"
#include "../common/types.hpp"
#include "../fiber/piecewise_linear_continuous_scalar_basis.hpp"

#include <map>
#include <memory>
#include <tbb/mutex.h>

namespace Bempp
{

/** \cond FORWARD_DECL */
class GridView;
/** \endcond */

/** \ingroup space
 *  \brief Space of continuous, piecewise linear scalar functions. */
template <typename BasisFunctionType>
class PiecewiseLinearScalarSpace : public ScalarSpace<BasisFunctionType>
{
public:
    typedef typename Space<BasisFunctionType>::CoordinateType CoordinateType;
    typedef typename Space<BasisFunctionType>::ComplexType ComplexType;

    explicit PiecewiseLinearScalarSpace(const shared_ptr<const Grid>& grid);
    virtual ~PiecewiseLinearScalarSpace();

    virtual int domainDimension() const;
    virtual int codomainDimension() const;


    /** \brief Return the variant of element \p element.
     *
     *  Possible return values:
     *    - 2: one-dimensional segment,
     *    - 3: triangular element,
     *    - 4: quadrilateral element. */
    virtual ElementVariant elementVariant(const Entity<0>& element) const;
    virtual void setElementVariant(const Entity<0>& element,
                                   ElementVariant variant);

    virtual const Fiber::Basis<BasisFunctionType>& basis(const Entity<0>& element) const;

private:
    /** \cond PRIVATE */
    Fiber::PiecewiseLinearContinuousScalarBasis<2, BasisFunctionType> m_lineBasis;
    Fiber::PiecewiseLinearContinuousScalarBasis<3, BasisFunctionType> m_triangleBasis;
    Fiber::PiecewiseLinearContinuousScalarBasis<4, BasisFunctionType> m_quadrilateralBasis;
    /** \endcond */
};

} // namespace Bempp

#endif
