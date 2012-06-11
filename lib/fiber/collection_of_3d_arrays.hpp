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

#ifndef fiber_collection_of_3d_arrays_hpp
#define fiber_collection_of_3d_arrays_hpp

#include "_3d_array.hpp"

#include <boost/scoped_array.hpp>

namespace Fiber
{

template <typename T> class CollectionOf2dSlicesOf3dArrays;
template <typename T> class CollectionOf1dSlicesOf3dArrays;
template <typename T> class CollectionOf2dSlicesOfConst3dArrays;
template <typename T> class CollectionOf1dSlicesOfConst3dArrays;

template <typename T>
class CollectionOf3dArrays
{
public:
    CollectionOf3dArrays();
    explicit CollectionOf3dArrays(int arrayCount);

    void set_size(int new_size);
    int size() const;

    _3dArray<T>& operator[](int index);
    const _3dArray<T>& operator[](int index) const;

    void fill(const T& value);

    CollectionOf2dSlicesOf3dArrays<T> slice(int index2);
    CollectionOf2dSlicesOfConst3dArrays<T> const_slice(int index2) const;
    CollectionOf1dSlicesOf3dArrays<T> slice(int index1, int index2);
    CollectionOf1dSlicesOfConst3dArrays<T> const_slice(int index1, int index2) const;

private:
    _3dArray<T>& array(int index);
    const _3dArray<T>& array(int index) const;
    void check_array_index(int array_index) const;

private:
    // Disable copy constructor and assignment operator
    CollectionOf3dArrays(const CollectionOf3dArrays& rhs);
    CollectionOf3dArrays& operator=(const CollectionOf3dArrays& rhs);

private:
    int m_size;
    boost::scoped_array<_3dArray<T> > m_arrays;
};

template <typename T>
class CollectionOf2dSlicesOf3dArrays
{
public:
    CollectionOf2dSlicesOf3dArrays(CollectionOf3dArrays<T>& collection, int index2);

    /** \brief Returns a reference to self.

      Useful to make a temporary CollectionOf2dSlicesOf3dArrays<T> an rvalue
      and pass it to a function accepting a reference to a non-const
      CollectionOf2dSlicesOf3dArrays<T>.

      Once we switch to C++11, this function can be removed because of the new
      support for rvalue references. */
    CollectionOf2dSlicesOf3dArrays& self();

    _2dSliceOfConst3dArray<T> operator[](int index) const;
    _2dSliceOf3dArray<T> operator[](int index);

    int size() const;

private:
    CollectionOf3dArrays<T>& m_collection;
    int m_index2;
};

template <typename T>
class CollectionOf2dSlicesOfConst3dArrays
{
public:
    CollectionOf2dSlicesOfConst3dArrays(
            const CollectionOf3dArrays<T>& collection, int index2);

    _2dSliceOfConst3dArray<T> operator[](int index) const;

    int size() const;

private:
    const CollectionOf3dArrays<T>& m_collection;
    int m_index2;
};

template <typename T>
class CollectionOf1dSlicesOf3dArrays
{
public:
    CollectionOf1dSlicesOf3dArrays(CollectionOf3dArrays<T>& collection,
                                   int index1, int index2);

    /** \brief Returns a reference to self.

      Useful to make a temporary CollectionOf1dSlicesOf3dArrays<T> an rvalue
      and pass it to a function accepting a reference to a non-const
      CollectionOf1dSlicesOf3dArrays<T>.

      Once we switch to C++11, this function can be removed because of the new
      support for rvalue references. */
    CollectionOf1dSlicesOf3dArrays& self();

    _1dSliceOfConst3dArray<T> operator[](int index) const;
    _1dSliceOf3dArray<T> operator[](int index);

    int size() const;

private:
    CollectionOf3dArrays<T>& m_collection;
    int m_index1, m_index2;
};

template <typename T>
class CollectionOf1dSlicesOfConst3dArrays
{
public:
    CollectionOf1dSlicesOfConst3dArrays(
            const CollectionOf3dArrays<T>& collection, int index1, int index2);

    _1dSliceOfConst3dArray<T> operator[](int index) const;

    int size() const;

private:
    const CollectionOf3dArrays<T>& m_collection;
    int m_index1, m_index2;
};

} // namespace Fiber

#include "collection_of_3d_arrays_imp.hpp"

#endif
