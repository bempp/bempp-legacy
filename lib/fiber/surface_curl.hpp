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

#ifndef fiber__surface_curl_hpp
#define fiber__surface_curl_hpp

#include "expression.hpp"

#include <armadillo>
#include <vector>

namespace Fiber {

template <typename ValueType>
class SurfaceCurl : public Expression<ValueType>
{
public:
    SurfaceCurl(int gridDimension) :
        m_gridDimension(gridDimension)
    {}

    virtual int domainDimension() const {
        return 1;
    }

    virtual int codomainDimension() const {
        return m_gridDimension;
    }

    virtual void addDependencies(int& basisDeps, int& geomDeps) const {
        addShapeFunctionSurfaceCurlDependencies(basisDeps, geomDeps);
    }

    virtual void evaluate(const BasisData<ValueType>& basisData,
                          const GeometricalData<ValueType>& geomData,
                          arma::Cube<ValueType>& result) const {
        // Probably the compatibility of dimensions of geomData's elements with
        // m_gridDimension should be checked.
        ScalarSpaceMapping::evaluateShapeFunctionSurfaceCurlsInternal(
                    basisData, geomData, result);
    }

private:
    const int m_gridDimension;
};

// namespace Fiber

#endif
