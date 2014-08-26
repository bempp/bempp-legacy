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

#include "vector.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <iostream>

#ifdef WITH_TRILINOS
#include <Teuchos_ArrayRCP.hpp>
#endif

namespace Bempp {

template <typename ValueType>
Vector<ValueType>::Vector(const arma::Col<ValueType> &vec) {
#ifdef WITH_TRILINOS
  const size_t size = vec.n_rows;
  Teuchos::ArrayRCP<ValueType> data(size);
  for (size_t i = 0; i < size; ++i)
    data[i] = vec(i);
  this->initialize(Thyra::defaultSpmdVectorSpace<ValueType>(size), data,
                   1 /* stride */);
#else
  m_vec = vec;
#endif
}

template <typename ValueType> Vector<ValueType>::Vector(size_t n) {
#ifdef WITH_TRILINOS
  Teuchos::ArrayRCP<ValueType> data(n);
  this->initialize(Thyra::defaultSpmdVectorSpace<ValueType>(n), data,
                   1 /* stride */);
#else
  m_vec.set_size(n);
#endif
}

template <typename ValueType> size_t Vector<ValueType>::size() const {
#ifdef WITH_TRILINOS
  return this->range()->dim();
#else
  return m_vec.n_rows;
#endif
}

template <typename ValueType> void Vector<ValueType>::dump() const {
#ifdef WITH_TRILINOS
  std::cout << asArmadilloVector() << std::endl; // inefficient
#else
  std::cout << m_vec << std::endl;
#endif
}

template <typename ValueType>
arma::Col<ValueType> Vector<ValueType>::asArmadilloVector() const {
#ifdef WITH_TRILINOS
  const size_t size = this->range()->dim();
  arma::Col<ValueType> col(size);
  for (size_t i = 0; i < size; ++i)
    col(i) = this->getPtr()[i];
  return col;
#else
  return m_vec;
#endif
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(Vector);

} // namespace Bempp
