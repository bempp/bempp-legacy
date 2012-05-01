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

#ifndef fiber_surface_curl_3d_hpp
#define fiber_surface_curl_3d_hpp

#include "expression.hpp"
#include "scalar_space_mapping.hpp"

#include <armadillo>
#include <vector>

namespace Fiber {

template <typename CoordinateType>
class SurfaceCurl3D : public Expression<CoordinateType>
{
public:    
    typedef typename Expression<CoordinateType>::ComplexType ComplexType;

    virtual int domainDimension() const {
        return 1;
    }

    virtual int codomainDimension() const {
        return 3;
    }

    virtual void addDependencies(int& basisDeps, int& geomDeps) const {
        ScalarSpaceMapping<CoordinateType>::
                addSurfaceCurlDependencies(basisDeps, geomDeps);
    }

private:
    virtual void evaluateImplReal(const BasisData<CoordinateType>& basisData,
                                  const GeometricalData<CoordinateType>& geomData,
                                  arma::Cube<CoordinateType>& result) const {
        ScalarSpaceMapping<CoordinateType>::
                evaluateSurfaceCurls3D(basisData, geomData, result);
    }

    virtual void evaluateImplComplex(const BasisData<ComplexType>& basisData,
                                     const GeometricalData<CoordinateType>& geomData,
                                     arma::Cube<ComplexType>& result) const {
        ScalarSpaceMapping<ComplexType>::
                evaluateSurfaceCurls3D(basisData, geomData, result);
    }
};

} // namespace Fiber

#endif
