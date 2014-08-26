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

#include "../common/common.hpp"
#include "boost/operators.hpp"

#include <stdexcept>

#ifndef NDEBUG
#define FIBER_CHECK_ARRAY_BOUNDS
#endif

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename T> class _4dSliceOf4dArray;
template <typename T> class _2dSliceOf4dArray;
template <typename T> class _1dSliceOf4dArray;
template <typename T> class Const4dSliceOf4dArray;
template <typename T> class Const2dSliceOf4dArray;
template <typename T> class Const1dSliceOf4dArray;
/** \endcond */

/** \brief Simple implementation of a 4D Fortran-ordered array.

Bound checking can optionally be activated by defining the symbol
FIBER_CHECK_ARRAY_BOUNDS. */
template <typename T>
class _4dArray
    : boost::additive<_4dArray<T>, boost::multiplicative<_4dArray<T>, T>> {
public:
  _4dArray();
  _4dArray(size_t extent0, size_t extent1, size_t extent2, size_t extent3);
  _4dArray(size_t extent0, size_t extent1, size_t extent2, size_t extent3,
           T *data);
  _4dArray(const _4dArray &rhs);
  _4dArray &operator=(const _4dArray &rhs);

  ~_4dArray();

  T &operator()(size_t index0, size_t index1, size_t index2, size_t index3);
  const T &operator()(size_t index0, size_t index1, size_t index2,
                      size_t index3) const;

  size_t extent(size_t dimension) const;
  void set_size(size_t extent0, size_t extent1, size_t extent2, size_t extent3);

  _4dArray<T> &operator+=(const _4dArray<T> &other);
  _4dArray<T> &operator*=(const T &other);

  typedef T *iterator;
  typedef const T *const_iterator;

  iterator begin();
  const_iterator begin() const;
  iterator end();
  const_iterator end() const;

private:
  void init_memory(size_t extent0, size_t extent1, size_t extent2,
                   size_t extent3);
  void free_memory();

#ifdef FIBER_CHECK_ARRAY_BOUNDS
  void check_dimension(size_t dimension) const;
  void check_extents(size_t extent0, size_t extent1, size_t extent2,
                     size_t extent3) const;
  void check_indices(size_t index0, size_t index1, size_t index2,
                     size_t index3) const;
#endif

private:
  size_t m_extents[4];
  bool m_owns;
  T *m_storage;
};

/** \brief Lightweight encapsulation of a 3D slice of a 4D array. */
template <typename T> class _3dSliceOf4dArray {
public:
  _3dSliceOf4dArray(_4dArray<T> &array, size_t index3);

  /** \brief Returns a reference to self.

    Useful to make a temporary _3dSliceOf4dArray<T> an rvalue and pass it to
    a function accepting a reference to a non-const _3dSliceOf4dArray<T>.

    Once we switch to C++11, this function can be removed because of the new
    support for rvalue references. */
  _3dSliceOf4dArray &self();

  const T &operator()(size_t index0, size_t index1, size_t index2) const;
  T &operator()(size_t index0, size_t index1, size_t index2);

  size_t extent(size_t dimension) const;

private:
  void check_dimension(size_t dimension) const;

private:
  _4dArray<T> &m_array;
  size_t m_index3;
};

/** \brief Lightweight encapsulation of a 2D slice of a constant 4d array. */
template <typename T> class _3dSliceOfConst4dArray {
public:
  _3dSliceOfConst4dArray(const _4dArray<T> &array, size_t index3);

  const T &operator()(size_t index0, size_t index1, size_t index2) const;

  size_t extent(size_t dimension) const;

private:
  void check_dimension(size_t dimension) const;

private:
  const _4dArray<T> &m_array;
  size_t m_index3;
};

/** \brief Lightweight encapsulation of a 2D slice of a 4d array. */
template <typename T> class _2dSliceOf4dArray {
public:
  _2dSliceOf4dArray(_4dArray<T> &array, size_t index2, size_t index3);

  /** \brief Returns a reference to self.

    Useful to make a temporary _2dSliceOf4dArray<T> an rvalue and pass it to
    a function accepting a reference to a non-const _2dSliceOf4dArray<T>.

    Once we switch to C++11, this function can be removed because of the new
    support for rvalue references. */
  _2dSliceOf4dArray &self();

  const T &operator()(size_t index0, size_t index1) const;
  T &operator()(size_t index0, size_t index1);

  size_t extent(size_t dimension) const;

private:
  void check_dimension(size_t dimension) const;

private:
  _4dArray<T> &m_array;
  size_t m_index2, m_index3;
};

/** \brief Lightweight encapsulation of a 2D slice of a constant 4d array. */
template <typename T> class _2dSliceOfConst4dArray {
public:
  _2dSliceOfConst4dArray(const _4dArray<T> &array, size_t index2,
                         size_t index3);

  const T &operator()(size_t index0, size_t index1) const;

  size_t extent(size_t dimension) const;

private:
  void check_dimension(size_t dimension) const;

private:
  const _4dArray<T> &m_array;
  size_t m_index2, m_index3;
};

/** \brief Lightweight encapsulation of a 1D slice of a 4d array. */
template <typename T> class _1dSliceOf4dArray {
public:
  /** \brief Construct a slice consisting of the elements
   * array(:,index1,index2,index3) */
  _1dSliceOf4dArray(_4dArray<T> &array, size_t index1, size_t index2,
                    size_t index3);

  /** \brief Returns a reference to self.

    Useful to make a temporary _1dSliceOf4dArray<T> an rvalue and pass it to
    a function accepting a reference to a non-const _2dSliceOf4dArray<T>.

    Once we switch to C++11, this function can be removed because of the new
    support for rvalue references. */
  _1dSliceOf4dArray &self();

  const T &operator()(size_t index0) const;
  T &operator()(size_t index0);

  size_t extent(size_t dimension) const;

private:
  void check_dimension(size_t dimension) const;

private:
  _4dArray<T> &m_array;
  size_t m_index1, m_index2, m_index3;
};

/** \brief Lightweight encapsulation of a 2D slice of a constant 4d array. */
template <typename T> class _1dSliceOfConst4dArray {
public:
  _1dSliceOfConst4dArray(const _4dArray<T> &array, size_t index1, size_t index2,
                         size_t index3);

  const T &operator()(size_t index0) const;

  size_t extent(size_t dimension) const;

private:
  void check_dimension(size_t dimension) const;

private:
  const _4dArray<T> &m_array;
  size_t m_index1, m_index2, m_index3;
};

} // namespace Fiber

#include "_4d_array_imp.hpp"

#endif
