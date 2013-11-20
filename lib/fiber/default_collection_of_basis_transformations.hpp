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

#ifndef fiber_default_collection_of_basis_transformations_hpp
#define fiber_default_collection_of_basis_transformations_hpp

#include "default_collection_of_shapeset_transformations.hpp"

namespace Fiber
{

/** \ingroup weak_form_elements
 *  \brief Default implementation of a collection of shape function transformations.
 *
 *  \deprecated This class is deprecated and provided only for compatibility
 *  reasons. Use DefaultCollectionOfShapesetTransformations instead.
 */
template <typename Functor>
class DefaultCollectionOfBasisTransformations :
        public DefaultCollectionOfShapesetTransformations<Functor>
{
public:
    explicit DefaultCollectionOfBasisTransformations(const Functor& functor) :
        DefaultCollectionOfShapesetTransformations<Functor>(functor)
    {}
};

} // namespace Fiber

#endif
