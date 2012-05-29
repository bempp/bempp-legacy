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

#ifndef fiber_scalar_function_value_hpp
#define fiber_scalar_function_value_hpp

#include "expression.hpp"
#include "scalar_space_mapping.hpp"
#include "CL/scalar_function_value.cl.str"

#include <armadillo>

namespace Fiber
{

template <typename CoordinateType>
class ScalarFunctionValue : public Expression<CoordinateType>
{
public:
    typedef typename Expression<CoordinateType>::ComplexType ComplexType;

    virtual int domainDimension() const {
        return 1;
    }

    virtual int codomainDimension() const {
        return 1;
    }

    virtual void addDependencies(int& basisDeps, int& geomDeps) const {
        ScalarSpaceMapping<CoordinateType>::
                addShapeFunctionDependencies(basisDeps, geomDeps);
    }

    // This seems unused
    virtual const std::string clStringEvaluate (const std::string modifier)
        const {
        std::string funcName ("devExpressionEvaluate");
        std::string str (scalar_function_value_cl,
             scalar_function_value_cl_len);
        if (modifier.size() > 0) {
            int n = str.find (funcName);
            if (n != std::string::npos) {
                size_t len = funcName.size();
                funcName.append (modifier);
                str.replace (n, len, funcName);
            }
        }
        return str;
    }

private:
    virtual void evaluateImplReal(const BasisData<CoordinateType>& basisData,
                                  const GeometricalData<CoordinateType>& geomData,
                                  arma::Cube<CoordinateType>& result) const {
        ScalarSpaceMapping<CoordinateType>::
                evaluateShapeFunctions(basisData, geomData, result);
    }

    virtual void evaluateImplComplex(const BasisData<ComplexType>& basisData,
                                     const GeometricalData<CoordinateType>& geomData,
                                     arma::Cube<ComplexType>& result) const {
        ScalarSpaceMapping<ComplexType>::
                evaluateShapeFunctions(basisData, geomData, result);
    }
};

} // namespace Fiber

#endif
