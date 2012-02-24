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

#ifndef bempp_linear_operator_superposition_hpp
#define bempp_linear_operator_superposition_hpp

#include "linear_operator.hpp"

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/tuple/tuple.hpp>

namespace Bempp
{

// only scalar multipliers allowed, tensor ones would
// require knowledge of vector components distribution
// in the discrete operator
template <typename ValueType>
class LinearOperatorSuperposition : public LinearOperator<ValueType>
{
public:
    typedef typename LinearOperator<ValueType>::IntegrationManagerFactory
    IntegrationManagerFactory;
    typedef typename LinearOperator<ValueType>::IntegrationManager IntegrationManager;

    /* acquires ownership of the operators passed via terms */
    LinearOperatorSuperposition(
            boost::ptr_vector<LinearOperator<ValueType> >& terms,
            const std::vector<ValueType>& multipliers);

    // acquires ownership of the operators passed via terms
    LinearOperatorSuperposition(
            const boost::tuple<LinearOperator<ValueType>*,
            LinearOperator<ValueType>*>& terms,
            const boost::tuple<ValueType, ValueType>& multipliers);
    // possibly add variants for longer parameter lists

    virtual int trialComponentCount() const;
    virtual int testComponentCount() const;

    virtual std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
    assembleOperator(
            const arma::Mat<ctype>& testPoints,
            const Space<ValueType>& trialSpace,
            const IntegrationManagerFactory& factory,
            const AssemblyOptions& options) const;

    virtual std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
    assembleWeakForm(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            const IntegrationManagerFactory& factory,
            const AssemblyOptions& options) const;

private:
    void init(boost::ptr_vector<LinearOperator<ValueType> >& terms,
              const std::vector<ValueType>& multipliers);

private:
    boost::ptr_vector<LinearOperator<ValueType> > m_terms;
    std::vector<ValueType> m_multipliers;
};

} //namespace Bempp

#endif
