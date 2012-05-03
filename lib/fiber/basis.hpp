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

#ifndef fiber_basis_hpp
#define fiber_basis_hpp

#include "scalar_traits.hpp"
#include "types.hpp"

#include <armadillo>

namespace Fiber
{

template <typename ValueType> struct BasisData;

template <typename ValueType>
class Basis
{
public:
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    virtual ~Basis() {}

    virtual int size() const = 0;
    /** \brief Maximum polynomial order of basis elements. */
    virtual int order() const = 0;
    virtual void evaluate(int what,
                          const arma::Mat<CoordinateType>& points,
                          LocalDofIndex localDofIndex,
                          BasisData<ValueType>& data) const = 0;

    /**
     * \brief Returns an OpenCL code snippet for basis function evaluation
     * \note The code snippet must provide device function devBasisEval
     */
    virtual std::pair<const char*,int> clCodeString (bool isTestBasis) const {
        throw std::runtime_error("Basis: clCodeString not implemented yet");
    }
};

} // namespace Fiber

#endif
