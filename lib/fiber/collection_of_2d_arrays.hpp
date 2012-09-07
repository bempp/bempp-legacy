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

#ifndef fiber_collection_of_2d_arrays_hpp
#define fiber_collection_of_2d_arrays_hpp

#include "_2d_array.hpp"

#include <boost/scoped_array.hpp>

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename T> class CollectionOf1dSlicesOf2dArrays;
template <typename T> class CollectionOf1dSlicesOfConst2dArrays;
/** \endcond */

template <typename T>
class CollectionOf2dArrays
{
public:
    CollectionOf2dArrays();
    explicit CollectionOf2dArrays(size_t arrayCount);
    CollectionOf2dArrays(const CollectionOf2dArrays& rhs);
    CollectionOf2dArrays& operator=(const CollectionOf2dArrays& rhs);

    void set_size(size_t new_size);
    size_t size() const;

    _2dArray<T>& operator[](size_t index);
    const _2dArray<T>& operator[](size_t index) const;

    void fill(const T& value);

    CollectionOf1dSlicesOf2dArrays<T> slice(size_t index1);
    CollectionOf1dSlicesOfConst2dArrays<T> const_slice(size_t index1) const;

private:
    _2dArray<T>& array(size_t index);
    const _2dArray<T>& array(size_t index) const;
    void check_array_index(size_t array_index) const;

private:
    size_t m_size;
    boost::scoped_array<_2dArray<T> > m_arrays;
};


template <typename T>
class CollectionOf1dSlicesOf2dArrays
{
public:
    CollectionOf1dSlicesOf2dArrays(CollectionOf2dArrays<T>& collection,
                                   size_t index1);

    /** \brief Returns a reference to self.

      Useful to make a temporary CollectionOf1dSlicesOf2dArrays<T> an rvalue
      and pass it to a function accepting a reference to a non-const
      CollectionOf1dSlicesOf2dArrays<T>.

      Once we switch to C++11, this function can be removed because of the new
      support for rvalue references. */
    CollectionOf1dSlicesOf2dArrays& self();

    _1dSliceOfConst2dArray<T> operator[](size_t index) const;
    _1dSliceOf2dArray<T> operator[](size_t index);

    size_t size() const;

private:
    CollectionOf2dArrays<T>& m_collection;
    size_t m_index1;
};

template <typename T>
class CollectionOf1dSlicesOfConst2dArrays
{
public:
    CollectionOf1dSlicesOfConst2dArrays(
            const CollectionOf2dArrays<T>& collection, size_t index1);

    _1dSliceOfConst2dArray<T> operator[](size_t index) const;

    size_t size() const;

private:
    const CollectionOf2dArrays<T>& m_collection;
    size_t m_index1;
};

} // namespace Fiber

#include "collection_of_2d_arrays_imp.hpp"

#endif
