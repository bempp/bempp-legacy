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

#ifndef bempp_acc_hpp
#define bempp_acc_hpp

/** \file acc.hpp Bounds checking for STL containers. */

#include "common.hpp"

namespace Bempp
{

/** \brief Returns \p i'th element of the container \p v, optionally checking 
    bounds.

    If \c NDEBUG is not defined, this function returns <tt>v.at(i)</tt>,
    otherwise <tt>v[i]</tt>.
*/
template <typename Container>
inline typename Container::reference acc(Container& v, typename Container::size_type i)
{
#ifdef NDEBUG
    return v[i];
#else
    return v.at(i);
#endif
}

/** \brief Returns \p i'th element of the container \p v, optionally checking 
    bounds.

    If \c NDEBUG is not defined, this function returns <tt>v.at(i)</tt>,
    otherwise <tt>v[i]</tt>.
*/
template <typename Container>
inline typename Container::const_reference acc(const Container& v, typename Container::size_type i)
{
#ifdef NDEBUG
    return v[i];
#else
    return v.at(i);
#endif
}

} // namespace Bempp

#endif
