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

#ifndef bempp_scalar_space_hpp
#define bempp_scalar_space_hpp

#include "../common/common.hpp"

#include "space.hpp"

#include <boost/scoped_ptr.hpp>

namespace Bempp
{

/** \ingroup space
 *  \brief Base class for spaces of scalar-valued functions. */
template <typename BasisFunctionType>
class ScalarSpace : public Space<BasisFunctionType>
{
    typedef Space<BasisFunctionType> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::CollectionOfShapesetTransformations
    CollectionOfShapesetTransformations;
    typedef typename Base::CollectionOfBasisTransformations
    CollectionOfBasisTransformations;

    explicit ScalarSpace(const shared_ptr<const Grid>& grid);
    /** \brief Copy constructor. */
    ScalarSpace(const ScalarSpace& other);
    virtual ~ScalarSpace();
    /** \brief Assignment operator. */
    ScalarSpace& operator=(const ScalarSpace& other);

    virtual const CollectionOfShapesetTransformations&
    basisFunctionValue() const;

private:
    struct Impl;
    boost::scoped_ptr<Impl> m_impl;
};

} // namespace Bempp

#endif
