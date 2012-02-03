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

#ifndef bempp_function_family_adapter_hpp
#define bempp_function_family_adapter_hpp

#include "function_family.hpp"
#include "../fiber/function_family.hpp"
#include "../grid/common.hpp"
#include <armadillo>

namespace Bempp {

/** \brief Wrapper of Bempp::FunctionFamily conforming to the
  Fiber::FunctionFamily interface. */
template <typename ValueType>
class FunctionFamilyAdapter :
        public ::Fiber::FunctionFamily<ValueType, ctype,
        FunctionFamilyAdapter<ValueType> >
{
public:
    FunctionFamilyAdapter(FunctionFamily<ValueType>& family) :
        m_wrapped(family) {
    }

    int domainDimension() const {
        return asImp().domainDimension();
    }

     int worldDimension() const {
        return asImp().domainDimension() + 1;
    }

    int codomainDimension() const {
        return asImp().codomainDimension();
    }

    int size() const {
        return asImp().size();
    }

    bool needsJacobianInverseTransposed() const {
        return asImp().needsJacobianInverseTransposed();
    }

    void setEvaluationPoints(int pointCount,
                             const ctype* evaluationPoints) {
        const arma::Mat<ctype> arma_evaluationPoints(
                    const_cast<ctype*>(evaluationPoints),
                    domainDimension(), pointCount, false /* copy_mem */);
        asImp().setEvaluationPoints(arma_evaluationPoints);
    }

    int evaluationPointCount() const {
        return asImp().evaluationPointCount();
    }

    void evaluate(const ctype* jacobianInversesTransposed,
                  ValueType* result) const {
        const int pointCount = evaluationPointCount();
        const arma::Cube<ctype> arma_jacobianInversesTransposed(
                    const_cast<ctype*>(jacobianInversesTransposed),
                    worldDimension(), domainDimension(), pointCount,
                    false /* copy_mem */);
        arma::Cube<ValueType> arma_result(
                    result, codomainDimension(), size(), pointCount, false);
        asImp().evaluate(arma_jacobianInversesTransposed, arma_result);
    }

private:
    FunctionFamily<ValueType>& asImp() {
        return m_wrapped;
    }

    const FunctionFamily<ValueType>& asImp() const {
        return m_wrapped;
    }

    FunctionFamily<ValueType>& m_wrapped;
};

} // namespace Bempp

#endif
