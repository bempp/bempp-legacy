// Copyright (C) 2011-2012 by the Fiber Authors
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

#ifndef fiber_numerical_double_integrator_hpp
#define fiber_numerical_double_integrator_hpp

#include "double_integrator.hpp"

namespace Fiber
{

template <typename ValueType> class Expression;
template <typename ValueType> class Kernel;

/** \brief Integration over pairs of elements. */
template <typename ValueType, typename GeometryFactory>
class NumericalDoubleIntegrator :
        public DoubleIntegrator<ValueType>
{
public:
    typedef typename DoubleIntegrator<ValueType>::ElementIndexPair
    ElementIndexPair;

protected:
    static void setupGeometry(int elementIndex,
                              const arma::Mat<ValueType>& vertices,
                              const arma::Mat<int>& elementCornerIndices,
                              const arma::Mat<char>& auxElementData,
                              typename GeometryFactory::Geometry& geometry);
};

} // namespace Fiber

#include "numerical_double_integrator_imp.hpp"

#endif
