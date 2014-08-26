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

#include "collection_of_4d_arrays.hpp"

namespace Fiber {

// CollectionOf4dArrays

template <typename T>
inline CollectionOf4dArrays<T>::CollectionOf4dArrays()
    : m_size(0) {}

template <typename T>
inline CollectionOf4dArrays<T>::CollectionOf4dArrays(size_t size)
    : m_size(size), m_arrays(new _4dArray<T>[size]) {
  // should we initialise to 0?
}

template <typename T>
inline void CollectionOf4dArrays<T>::set_size(size_t new_size) {
  if (new_size == m_size)
    return;
  m_arrays.reset(new_size == 0 ? 0 : new _4dArray<T>[new_size]);
  m_size = new_size;
}

template <typename T> inline size_t CollectionOf4dArrays<T>::size() const {
  return m_size;
}

template <typename T>
inline void CollectionOf4dArrays<T>::fill(const T &value) {
  for (size_t a = 0; a < m_size; ++a)
    std::fill((*this)[a].begin(), (*this)[a].end(), value);
}

template <typename T>
inline _4dArray<T> &CollectionOf4dArrays<T>::operator[](size_t index) {
  return array(index);
}

template <typename T>
inline const _4dArray<T> &CollectionOf4dArrays<T>::
operator[](size_t index) const {
  return array(index);
}

template <typename T>
inline CollectionOf3dSlicesOf4dArrays<T>
CollectionOf4dArrays<T>::slice(size_t index3) {
  return CollectionOf3dSlicesOf4dArrays<T>(*this, index3);
}

template <typename T>
inline CollectionOf3dSlicesOfConst4dArrays<T>
CollectionOf4dArrays<T>::const_slice(size_t index3) const {
  return CollectionOf3dSlicesOfConst4dArrays<T>(*this, index3);
}

template <typename T>
inline CollectionOf2dSlicesOf4dArrays<T>
CollectionOf4dArrays<T>::slice(size_t index2, size_t index3) {
  return CollectionOf2dSlicesOf4dArrays<T>(*this, index2, index3);
}

template <typename T>
inline CollectionOf2dSlicesOfConst4dArrays<T>
CollectionOf4dArrays<T>::const_slice(size_t index2, size_t index3) const {
  return CollectionOf2dSlicesOfConst4dArrays<T>(*this, index2, index3);
}

template <typename T>
inline CollectionOf1dSlicesOf4dArrays<T>
CollectionOf4dArrays<T>::slice(size_t index1, size_t index2, size_t index3) {
  return CollectionOf1dSlicesOf4dArrays<T>(*this, index1, index2, index3);
}

template <typename T>
inline CollectionOf1dSlicesOfConst4dArrays<T>
CollectionOf4dArrays<T>::const_slice(size_t index1, size_t index2,
                                     size_t index3) const {
  return CollectionOf1dSlicesOfConst4dArrays<T>(*this, index1, index2, index3);
}

template <typename T>
inline _4dArray<T> &CollectionOf4dArrays<T>::array(size_t i) {
#ifdef FIBER_CHECK_ARRAY_BOUNDS
  check_array_index(i);
#endif
  return m_arrays[i];
}

template <typename T>
inline const _4dArray<T> &CollectionOf4dArrays<T>::array(size_t i) const {
#ifdef FIBER_CHECK_ARRAY_BOUNDS
  check_array_index(i);
#endif
  return m_arrays[i];
}

template <typename T>
inline void
CollectionOf4dArrays<T>::check_array_index(size_t array_index) const {
#ifdef FIBER_CHECK_ARRAY_BOUNDS
  if (m_size <= array_index)
    throw std::invalid_argument("Invalid array index");
#endif
}

// CollectionOf3dSlicesOf4dArrays

template <typename T>
inline CollectionOf3dSlicesOf4dArrays<T>::CollectionOf3dSlicesOf4dArrays(
    CollectionOf4dArrays<T> &collection, size_t index3)
    : m_collection(collection), m_index3(index3) {}

template <typename T>
inline CollectionOf3dSlicesOf4dArrays<T> &
CollectionOf3dSlicesOf4dArrays<T>::self() {
  return *this;
}

template <typename T>
inline _3dSliceOfConst4dArray<T> CollectionOf3dSlicesOf4dArrays<T>::
operator[](size_t index) const {
  return _3dSliceOfConst4dArray<T>(m_collection[index], m_index3);
}

template <typename T>
inline _3dSliceOf4dArray<T> CollectionOf3dSlicesOf4dArrays<T>::
operator[](size_t index) {
  return _3dSliceOf4dArray<T>(m_collection[index], m_index3);
}

template <typename T>
inline size_t CollectionOf3dSlicesOf4dArrays<T>::size() const {
  return m_collection.size();
}

// CollectionOf3dSlicesOfConst4dArrays

template <typename T>
inline CollectionOf3dSlicesOfConst4dArrays<
    T>::CollectionOf3dSlicesOfConst4dArrays(const CollectionOf4dArrays<T> &
                                                collection,
                                            size_t index3)
    : m_collection(collection), m_index3(index3) {}

template <typename T>
inline _3dSliceOfConst4dArray<T> CollectionOf3dSlicesOfConst4dArrays<T>::
operator[](size_t index) const {
  return _3dSliceOfConst4dArray<T>(m_collection[index], m_index3);
}

template <typename T>
inline size_t CollectionOf3dSlicesOfConst4dArrays<T>::size() const {
  return m_collection.size();
}

// CollectionOf2dSlicesOf4dArrays

template <typename T>
inline CollectionOf2dSlicesOf4dArrays<T>::CollectionOf2dSlicesOf4dArrays(
    CollectionOf4dArrays<T> &collection, size_t index2, size_t index3)
    : m_collection(collection), m_index2(index2), m_index3(index3) {}

template <typename T>
inline CollectionOf2dSlicesOf4dArrays<T> &
CollectionOf2dSlicesOf4dArrays<T>::self() {
  return *this;
}

template <typename T>
inline _2dSliceOfConst4dArray<T> CollectionOf2dSlicesOf4dArrays<T>::
operator[](size_t index) const {
  return _2dSliceOfConst4dArray<T>(m_collection[index], m_index2, m_index3);
}

template <typename T>
inline _2dSliceOf4dArray<T> CollectionOf2dSlicesOf4dArrays<T>::
operator[](size_t index) {
  return _2dSliceOf4dArray<T>(m_collection[index], m_index2, m_index3);
}

template <typename T>
inline size_t CollectionOf2dSlicesOf4dArrays<T>::size() const {
  return m_collection.size();
}

// CollectionOf2dSlicesOfConst4dArrays

template <typename T>
inline CollectionOf2dSlicesOfConst4dArrays<
    T>::CollectionOf2dSlicesOfConst4dArrays(const CollectionOf4dArrays<T> &
                                                collection,
                                            size_t index2, size_t index3)
    : m_collection(collection), m_index2(index2), m_index3(index3) {}

template <typename T>
inline _2dSliceOfConst4dArray<T> CollectionOf2dSlicesOfConst4dArrays<T>::
operator[](size_t index) const {
  return _2dSliceOfConst4dArray<T>(m_collection[index], m_index2, m_index3);
}

template <typename T>
inline size_t CollectionOf2dSlicesOfConst4dArrays<T>::size() const {
  return m_collection.size();
}

// CollectionOf1dSlicesOf4dArrays

template <typename T>
inline CollectionOf1dSlicesOf4dArrays<T>::CollectionOf1dSlicesOf4dArrays(
    CollectionOf4dArrays<T> &collection, size_t index1, size_t index2,
    size_t index3)
    : m_collection(collection), m_index1(index1), m_index2(index2),
      m_index3(index3) {}

template <typename T>
inline CollectionOf1dSlicesOf4dArrays<T> &
CollectionOf1dSlicesOf4dArrays<T>::self() {
  return *this;
}

template <typename T>
inline _1dSliceOfConst4dArray<T> CollectionOf1dSlicesOf4dArrays<T>::
operator[](size_t index) const {
  return _1dSliceOfConst4dArray<T>(m_collection[index], m_index1, m_index2,
                                   m_index3);
}

template <typename T>
inline _1dSliceOf4dArray<T> CollectionOf1dSlicesOf4dArrays<T>::
operator[](size_t index) {
  return _1dSliceOf4dArray<T>(m_collection[index], m_index1, m_index2,
                              m_index3);
}

template <typename T>
inline size_t CollectionOf1dSlicesOf4dArrays<T>::size() const {
  return m_collection.size();
}

// CollectionOf1dSlicesOfConst4dArrays

template <typename T>
inline CollectionOf1dSlicesOfConst4dArrays<
    T>::CollectionOf1dSlicesOfConst4dArrays(const CollectionOf4dArrays<T> &
                                                collection,
                                            size_t index1, size_t index2,
                                            size_t index3)
    : m_collection(collection), m_index1(index1), m_index2(index2),
      m_index3(index3) {}

template <typename T>
inline _1dSliceOfConst4dArray<T> CollectionOf1dSlicesOfConst4dArrays<T>::
operator[](size_t index) const {
  return _1dSliceOfConst4dArray<T>(m_collection[index], m_index1, m_index2,
                                   m_index3);
}

template <typename T>
inline size_t CollectionOf1dSlicesOfConst4dArrays<T>::size() const {
  return m_collection.size();
}

} // namespace Fiber
