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

#ifndef bempp_function_family_hpp
#define bempp_function_family_hpp

#include <armadillo>

namespace Bempp {

template <typename ValueType>
class FunctionFamily
{
public:
    virtual ~FunctionFamily() {}

    virtual int domainDimension() const = 0;
    virtual int codomainDimension() const = 0;
    virtual bool needsJacobianInverseTransposed() const = 0;
    virtual int size() const = 0;

    virtual int evaluationPointCount() const = 0;

    virtual void setEvaluationPoints(const arma::Mat<ctype>& evaluationPoints) = 0;
    virtual void evaluate(const arma::Cube<ctype>& jacobianInverseTransposed,
                          arma::Cube<ValueType>& result) const = 0;
};

} // namespace Bempp

#endif
