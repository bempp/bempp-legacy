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

#ifndef bempp_vector_hpp
#define bempp_vector_hpp

#include "config_trilinos.hpp"

#include <armadillo>

#ifdef WITH_TRILINOS
#include <Thyra_DefaultSpmdVector.hpp>
#endif

namespace Bempp {

/** \brief Encapsulation of a vector.

  If BEM++ is compiled with Trilinos, Vector is implemented by means of
  Thyra::DefaultSpmdVector; otherwise an Armadillo-based fallback
  implementation is used. */
template <typename ValueType>
class Vector
#ifdef WITH_TRILINOS
        : public Thyra::DefaultSpmdVector<ValueType>
#endif
{
public:
    /** Construct the discrete source term from an Armadillo vector. */
    Vector(const arma::Col<ValueType>& vec);

    /** \brief Write a textual representation of the source term to standard output. */
    void dump() const;

    /** \brief Vector representation of the source term. */
    arma::Col<ValueType> asArmadilloVector() const;

private:
#ifndef WITH_TRILINOS
    // fallback implementation
    arma::Col<ValueType> m_vec;
#endif
};

} // namespace Bempp

#endif
