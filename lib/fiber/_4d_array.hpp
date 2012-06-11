// Copyright (C) 2011-2012 by the Bem++ Authors
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

#ifndef fiber_4d_array_hpp
#define fiber_4d_array_hpp

#include <stdexcept>

namespace Fiber {

template <typename T> class _4dSliceOf4dArray;
template <typename T> class _2dSliceOf4dArray;
template <typename T> class _1dSliceOf4dArray;
template <typename T> class Const4dSliceOf4dArray;
template <typename T> class Const2dSliceOf4dArray;
template <typename T> class Const1dSliceOf4dArray;

/** \brief Simple implementation of a 4D Fortran-ordered array.

Bound checking can optionally be activated by defining the symbol
FIBER_CHECK_ARRAY_BOUNDS. */
template <typename T>
class _4dArray
{
public:
    _4dArray();
    _4dArray(int extent0, int extent1, int extent2, int extent3);
    _4dArray(int extent0, int extent1, int extent2, int extent3, T* data);

    ~_4dArray();

    T& operator()(int index0, int index1, int index2, int index3);
    const T& operator()(int index0, int index1, int index2, int index3) const;

    int extent(int dimension) const;
    void set_size(int extent0, int extent1, int extent2, int extent3);

    typedef T* iterator;
    typedef const T* const_iterator;

    iterator begin();
    const_iterator begin() const;
    iterator end();
    const_iterator end() const;

private:
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    void check_dimension(int dimension) const;
    void check_extents(int extent0, int extent1, int extent2, int extent3) const;
    void check_indices(int index0, int index1, int index2, int index3) const;
#endif

private:
    // Disable copy constructor and assignment operator
    _4dArray(const _4dArray& rhs);
    _4dArray& operator=(const _4dArray& rhs);

private:
    int m_extents[4];
    bool m_owns;
    T* m_storage;
};

/** \brief Lightweight encapsulation of a 3D slice of a 4D array. */
template <typename T>
class _3dSliceOf4dArray
{
public:
    _3dSliceOf4dArray(_4dArray<T>& array, int index3);

    /** \brief Returns a reference to self.

      Useful to make a temporary _3dSliceOf4dArray<T> an rvalue and pass it to
      a function accepting a reference to a non-const _3dSliceOf4dArray<T>.

      Once we switch to C++11, this function can be removed because of the new
      support for rvalue references. */
    _3dSliceOf4dArray& self();

    const T& operator()(int index0, int index1, int index2) const;
    T& operator()(int index0, int index1, int index2);

    int extent(int dimension) const;

private:
    void check_dimension(int dimension) const;

private:
    _4dArray<T>& m_array;
    int m_index3;
};

/** \brief Lightweight encapsulation of a 2D slice of a constant 4d array. */
template <typename T>
class _3dSliceOfConst4dArray
{
public:
    _3dSliceOfConst4dArray(const _4dArray<T>& array, int index3);

    const T& operator()(int index0, int index1, int index2) const;

    int extent(int dimension) const;

private:
    void check_dimension(int dimension) const;

private:
    const _4dArray<T>& m_array;
    int m_index3;
};

/** \brief Lightweight encapsulation of a 2D slice of a 4d array. */
template <typename T>
class _2dSliceOf4dArray
{
public:
    _2dSliceOf4dArray(_4dArray<T>& array, int index2, int index3);

    /** \brief Returns a reference to self.

      Useful to make a temporary _2dSliceOf4dArray<T> an rvalue and pass it to
      a function accepting a reference to a non-const _2dSliceOf4dArray<T>.

      Once we switch to C++11, this function can be removed because of the new
      support for rvalue references. */
    _2dSliceOf4dArray& self();

    const T& operator()(int index0, int index1) const;
    T& operator()(int index0, int index1);

    int extent(int dimension) const;

private:
    void check_dimension(int dimension) const;

private:
    _4dArray<T>& m_array;
    int m_index2, m_index3;
};

/** \brief Lightweight encapsulation of a 2D slice of a constant 4d array. */
template <typename T>
class _2dSliceOfConst4dArray
{
public:
    _2dSliceOfConst4dArray(const _4dArray<T>& array, int index2, int index3);

    const T& operator()(int index0, int index1) const;

    int extent(int dimension) const;

private:
    void check_dimension(int dimension) const;

private:
    const _4dArray<T>& m_array;
    int m_index2, m_index3;
};

/** \brief Lightweight encapsulation of a 1D slice of a 4d array. */
template <typename T>
class _1dSliceOf4dArray
{
public:
    /** \brief Construct a slice consisting of the elements array(:,index1,index2,index3) */
    _1dSliceOf4dArray(_4dArray<T>& array, int index1, int index2, int index3);

    /** \brief Returns a reference to self.

      Useful to make a temporary _1dSliceOf4dArray<T> an rvalue and pass it to
      a function accepting a reference to a non-const _2dSliceOf4dArray<T>.

      Once we switch to C++11, this function can be removed because of the new
      support for rvalue references. */
    _1dSliceOf4dArray& self();

    const T& operator()(int index0) const;
    T& operator()(int index0);

    int extent(int dimension) const;

private:
    void check_dimension(int dimension) const;

private:
    _4dArray<T>& m_array;
    int m_index1, m_index2, m_index3;
};

/** \brief Lightweight encapsulation of a 2D slice of a constant 4d array. */
template <typename T>
class _1dSliceOfConst4dArray
{
public:
    _1dSliceOfConst4dArray(const _4dArray<T>& array, int index1, int index2, int index3);

    const T& operator()(int index0) const;

    int extent(int dimension) const;

private:
    void check_dimension(int dimension) const;

private:
    const _4dArray<T>& m_array;
    int m_index1, m_index2, m_index3;
};

} // namespace Fiber

#include "_4d_array_imp.hpp"

#endif
