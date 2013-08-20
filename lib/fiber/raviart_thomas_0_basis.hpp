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

#ifndef fiber_raviart_thomas_0_basis_hpp
#define fiber_raviart_thomas_0_basis_hpp

#include "raviart_thomas_0_shapeset.hpp"
#include "../common/deprecated.hpp"

namespace Fiber
{

/** \brief Shapeset composed of the lowest-order Raviart-Thomas functions.
 *
 *  \deprecated This class is deprecated, use RaviartThomas0Shapeset instead. */
template <int elementVertexCount, typename ValueType>
class BEMPP_DEPRECATED RaviartThomas0Basis :
        public RaviartThomas0Shapeset<elementVertexCount, ValueType>
{
};

} // namespace Fiber

#endif
