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

#ifndef fiber_numerical_quadrature_hpp
#define fiber_numerical_quadrature_hpp

#include "../common/common.hpp"

/** \file
 *
 *  Low-level functions filling arrays of quadrature points and weights. */

#include "double_quadrature_descriptor.hpp"
#include "single_quadrature_descriptor.hpp"
#include "types.hpp"

namespace Fiber {

/** \brief Retrieve points and weights for a quadrature over a single element.
 *
 *  \param[in] elementCornerCount
 *    Number of corners of the element to be integrated on.
 *  \param[in] accuracyOrder
 *    Accuracy order of the quadrature, i.e. its degree of exactness.
 *  \param[out] points
 *    Quadrature points.
 *  \param[out] weights
 *    Quadrature weights. */
template <typename ValueType>
void fillSingleQuadraturePointsAndWeights(int elementCornerCount,
                                          int accuracyOrder,
                                          Matrix<ValueType> &points,
                                          std::vector<ValueType> &weights);

template <typename ValueType>
void fillDoubleSingularQuadraturePointsAndWeights(
    const DoubleQuadratureDescriptor &desc, Matrix<ValueType> &testPoints,
    Matrix<ValueType> &trialPoints, std::vector<ValueType> &weights);

} // namespace Fiber

#endif
