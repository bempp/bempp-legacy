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

#ifndef bempp_quadrature_selector_hpp
#define bempp_quadrature_selector_hpp

#include "../common/multidimensional_arrays.hpp"
#include "../grid/common.hpp"
#include "../fiber/quadrature_selector.hpp"

#include <armadillo>

namespace Bempp
{

class GeometryAdapter;

template <typename ValueType>
class QuadratureSelector
{
public:    
    virtual void selectDoubleQuadratureRules(
            const std::vector<const GeometryAdapter*>& testGeometries,
            const std::vector<const GeometryAdapter*>& trialGeometries,
            Array2D<Fiber::QuadratureRule>& quadRules) const = 0;
    virtual void doubleQuadratureRulePointsAndWeights(
            Fiber::QuadratureRule quadRule,
            arma::Mat<ctype>& testQuadPoints,
            arma::Mat<ctype>& trialQuadPoints,
            std::vector<ValueType>& quadWeights) const = 0;
};

} // namespace Bempp

#endif
