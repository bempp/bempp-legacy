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

#ifndef fiber_collection_of_4d_arrays_hpp
#define fiber_collection_of_4d_arrays_hpp

#include "_4d_array.hpp"

#include <boost/scoped_array.hpp>

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename T> class CollectionOf3dSlicesOf4dArrays;
template <typename T> class CollectionOf2dSlicesOf4dArrays;
template <typename T> class CollectionOf1dSlicesOf4dArrays;
template <typename T> class CollectionOf3dSlicesOfConst4dArrays;
template <typename T> class CollectionOf2dSlicesOfConst4dArrays;
template <typename T> class CollectionOf1dSlicesOfConst4dArrays;
/** \endcond */

template <typename T>
class CollectionOf4dArrays
{
public:
    CollectionOf4dArrays();
    explicit CollectionOf4dArrays(size_t arrayCount);

    void set_size(size_t new_size);
    size_t size() const;

    _4dArray<T>& operator[](size_t index);
    const _4dArray<T>& operator[](size_t index) const;

    void fill(const T& value);

    CollectionOf3dSlicesOf4dArrays<T> slice(size_t index3);
    CollectionOf3dSlicesOfConst4dArrays<T> const_slice(size_t index3) const;
    CollectionOf2dSlicesOf4dArrays<T> slice(size_t index2, size_t index3);
    CollectionOf2dSlicesOfConst4dArrays<T> const_slice(size_t index2, size_t index3) const;
    CollectionOf1dSlicesOf4dArrays<T> slice(size_t index1, size_t index2, size_t index3);
    CollectionOf1dSlicesOfConst4dArrays<T> const_slice(size_t index1, size_t index2, size_t index3) const;

private:
    _4dArray<T>& array(size_t index);
    const _4dArray<T>& array(size_t index) const;
    void check_array_index(size_t array_index) const;

private:
    // Disable copy constructor and assignment operator
    CollectionOf4dArrays(const CollectionOf4dArrays& rhs);
    CollectionOf4dArrays& operator=(const CollectionOf4dArrays& rhs);

private:
    size_t m_size;
    boost::scoped_array<_4dArray<T> > m_arrays;
};

template <typename T>
class CollectionOf3dSlicesOf4dArrays
{
public:
    CollectionOf3dSlicesOf4dArrays(CollectionOf4dArrays<T>& collection, size_t index3);

    /** \brief Returns a reference to self.

      Useful to make a temporary CollectionOf3dSlicesOf4dArrays<T> an rvalue
      and pass it to a function accepting a reference to a non-const
      CollectionOf3dSlicesOf4dArrays<T>.

      Once we switch to C++11, this function can be removed because of the new
      support for rvalue references. */
    CollectionOf3dSlicesOf4dArrays& self();

    _3dSliceOfConst4dArray<T> operator[](size_t index) const;
    _3dSliceOf4dArray<T> operator[](size_t index);

    size_t size() const;

private:
    CollectionOf4dArrays<T>& m_collection;
    size_t m_index3;
};

template <typename T>
class CollectionOf3dSlicesOfConst4dArrays
{
public:
    CollectionOf3dSlicesOfConst4dArrays(
            const CollectionOf4dArrays<T>& collection, size_t index3);

    _3dSliceOfConst4dArray<T> operator[](size_t index) const;

    size_t size() const;

private:
    const CollectionOf4dArrays<T>& m_collection;
    size_t m_index3;
};

template <typename T>
class CollectionOf2dSlicesOf4dArrays
{
public:
    CollectionOf2dSlicesOf4dArrays(CollectionOf4dArrays<T>& collection, size_t index2, size_t index3);

    /** \brief Returns a reference to self.

      Useful to make a temporary CollectionOf2dSlicesOf4dArrays<T> an rvalue
      and pass it to a function accepting a reference to a non-const
      CollectionOf2dSlicesOf4dArrays<T>.

      Once we switch to C++11, this function can be removed because of the new
      support for rvalue references. */
    CollectionOf2dSlicesOf4dArrays& self();

    _2dSliceOfConst4dArray<T> operator[](size_t index) const;
    _2dSliceOf4dArray<T> operator[](size_t index);

    size_t size() const;

private:
    CollectionOf4dArrays<T>& m_collection;
    size_t m_index2, m_index3;
};

template <typename T>
class CollectionOf2dSlicesOfConst4dArrays
{
public:
    CollectionOf2dSlicesOfConst4dArrays(
            const CollectionOf4dArrays<T>& collection, size_t index2, size_t index3);

    _2dSliceOfConst4dArray<T> operator[](size_t index) const;

    size_t size() const;

private:
    const CollectionOf4dArrays<T>& m_collection;
    size_t m_index2, m_index3;
};

template <typename T>
class CollectionOf1dSlicesOf4dArrays
{
public:
    CollectionOf1dSlicesOf4dArrays(CollectionOf4dArrays<T>& collection,
                                   size_t index1, size_t index2, size_t index3);

    /** \brief Returns a reference to self.

      Useful to make a temporary CollectionOf1dSlicesOf4dArrays<T> an rvalue
      and pass it to a function accepting a reference to a non-const
      CollectionOf1dSlicesOf4dArrays<T>.

      Once we switch to C++11, this function can be removed because of the new
      support for rvalue references. */
    CollectionOf1dSlicesOf4dArrays& self();

    _1dSliceOfConst4dArray<T> operator[](size_t index) const;
    _1dSliceOf4dArray<T> operator[](size_t index);

    size_t size() const;

private:
    CollectionOf4dArrays<T>& m_collection;
    size_t m_index1, m_index2, m_index3;
};

template <typename T>
class CollectionOf1dSlicesOfConst4dArrays
{
public:
    CollectionOf1dSlicesOfConst4dArrays(
            const CollectionOf4dArrays<T>& collection, size_t index1, size_t index2, size_t index3);

    _1dSliceOfConst4dArray<T> operator[](size_t index) const;

    size_t size() const;

private:
    const CollectionOf4dArrays<T>& m_collection;
    size_t m_index1, m_index2, m_index3;
};

} // namespace Fiber

#include "collection_of_4d_arrays_imp.hpp"

#endif
