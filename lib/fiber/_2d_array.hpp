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

#ifndef fiber_2d_array_hpp
#define fiber_2d_array_hpp

#include "../common/common.hpp"

#include <stdexcept>

namespace Fiber {

/** \brief Simple implementation of a 2D Fortran-ordered array.

Bound checking can optionally be activated by defining the symbol
FIBER_CHECK_ARRAY_BOUNDS. */
template <typename T>
class _2dArray
{
public:
    _2dArray();
    _2dArray(size_t extent0, size_t extent1);
    _2dArray(size_t extent0, size_t extent1, T* data);
    _2dArray(const _2dArray& other);
    _2dArray& operator=(const _2dArray& rhs);

    ~_2dArray();

    T& operator()(size_t index0, size_t index1);
    const T& operator()(size_t index0, size_t index1) const;

    size_t extent(size_t dimension) const;
    void set_size(size_t extent0, size_t extent1);

    typedef T* iterator;
    typedef const T* const_iterator;

    iterator begin();
    const_iterator begin() const;
    iterator end();
    const_iterator end() const;

private:
    void init_memory(size_t extent0, size_t extent1);
    void free_memory();

#ifdef FIBER_CHECK_ARRAY_BOUNDS
    void check_dimension(size_t dimension) const;
    void check_extents(size_t extent0, size_t extent1) const;
    void check_indices(size_t index0, size_t index1) const;
#endif

private:
    size_t m_extents[2];
    bool m_owns;
    T* m_storage;
};

/** \brief Lightweight encapsulation of a 1D slice of a 2D array. */
template <typename T>
class _1dSliceOf2dArray
{
public:
    /** \brief Construct a slice consisting of the elements array(:,index1) */
    _1dSliceOf2dArray(_2dArray<T>& array, size_t index1);

    /** \brief Returns a reference to self.

      Useful to make a temporary _1dSliceOf2dArray<T> an rvalue and pass it to
      a function accepting a reference to a non-const _2dSliceOf2dArray<T>.

      Once we switch to C++11, this function can be removed because of the new
      support for rvalue references. */
    _1dSliceOf2dArray& self();

    const T& operator()(size_t index0) const;
    T& operator()(size_t index0);

    size_t extent(size_t dimension) const;

private:
    void check_dimension(size_t dimension) const;

private:
    _2dArray<T>& m_array;
    size_t m_index1;
};

/** \brief Lightweight encapsulation of a 1D slice of a constant 2D array. */
template <typename T>
class _1dSliceOfConst2dArray
{
public:
    _1dSliceOfConst2dArray(const _2dArray<T>& array, size_t index1);

    const T& operator()(size_t index0) const;

    size_t extent(size_t dimension) const;

private:
    void check_dimension(size_t dimension) const;

private:
    const _2dArray<T>& m_array;
    size_t m_index1;
};

} // namespace Fiber

#include "_2d_array_imp.hpp"

#endif
