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

#ifndef bempp_kernel_adapter_hpp
#define bempp_kernel_adapter_hpp

#include "kernel.hpp"
#include "../fiber/kernel.hpp"

namespace Bempp
{

/** \brief Wrapper of Bempp::Kernel conforming to the Fiber::Kernel interface. */
template <typename ValueType>
class KernelAdapter :
        public Fiber::Kernel<ValueType, ctype, KernelAdapter<ValueType> >
{
public:
    KernelAdapter(const Kernel<ValueType>& kernel) :
        m_wrapped(kernel) {
    }

    int worldDimension() const  {
        return asImp().worldDimension();
    }

    int codomainDimension() const  {
        return asImp().codomainDimension();
    }

    int domainDimension() const {
        return asImp().domainDimension();
    }

    bool needsTestNormal() const {
        return asImp().needsTestNormal();
    }

    bool needsTrialNormal() const {
        return asImp().needsTrialNormal();
    }

    void evaluate(int pointCount,
                  const ctype* testPoints,
                  const ctype* trialPoints,
                  const ctype* testNormals,
                  const ctype* trialNormals,
                  ValueType* result) const {
        const int worldDim = worldDimension();

        // The const_casts are there so that we can use the usual "const
        // arma::Mat<ctype>" type for constant matrices passed to
        // Kernel::evaluate() instead of the clumsier variant "const
        // arma::Mat<const ctype>"
        const arma::Mat<ctype> arma_testPoints(
                    const_cast<ctype*>(testPoints), worldDim, pointCount, false /* copy_mem */);
        const arma::Mat<ctype> arma_trialPoints(
                    const_cast<ctype*>(trialPoints), worldDim, pointCount, false);
        const arma::Mat<ctype> arma_testNormals(
                    const_cast<ctype*>(testNormals), worldDim, pointCount, false);
        const arma::Mat<ctype> arma_trialNormals(
                    const_cast<ctype*>(trialNormals), worldDim, pointCount, false);
        arma::Cube<ValueType> arma_result(
                    result, codomainDimension(), domainDimension(), pointCount, false);
        asImp().evaluate(arma_testPoints, arma_trialPoints,
                         arma_testNormals, arma_trialNormals, arma_result);
    }

private:
    const Kernel<ValueType>& asImp() const {
        return m_wrapped;
    }

    const Kernel<ValueType>& m_wrapped;
};

} // namespace Bempp

#endif
