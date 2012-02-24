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

#include "linear_operator_superposition.hpp"
#include "discrete_scalar_valued_linear_operator_superposition.hpp"
#include "discrete_vector_valued_linear_operator_superposition.hpp"

namespace Bempp
{

template <typename ValueType>
LinearOperatorSuperposition<ValueType>::LinearOperatorSuperposition(
        boost::ptr_vector<LinearOperator<ValueType> >& terms,
        const std::vector<ValueType>& multipliers)
{
   init(terms, multipliers);
}

template <typename ValueType>
LinearOperatorSuperposition<ValueType>::LinearOperatorSuperposition(
        const boost::tuple<LinearOperator<ValueType>*,
        LinearOperator<ValueType>*>& terms,
        const boost::tuple<ValueType, ValueType>& multipliers)
{
    boost::ptr_vector<LinearOperator<ValueType> > vTerms;
    vTerms.push_back(terms.template get<0>());
    vTerms.push_back(terms.template get<1>());
    std::vector<ValueType> vMultipliers;
    vMultipliers.push_back(multipliers.template get<0>());
    vMultipliers.push_back(multipliers.template get<1>());
    init(vTerms, vMultipliers);
}

template <typename ValueType>
int LinearOperatorSuperposition<ValueType>::trialComponentCount() const
{
    if (m_terms.empty())
        return 1;
    else
        return m_terms[0].trialComponentCount();
}

template <typename ValueType>
int LinearOperatorSuperposition<ValueType>::testComponentCount() const
{
    if (m_terms.empty())
        return 1;
    else
        return m_terms[0].testComponentCount();
}

template <typename ValueType>
std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
LinearOperatorSuperposition<ValueType>::assembleOperator(
        const arma::Mat<ctype>& testPoints,
        const Space<ValueType>& trialSpace,
        const IntegrationManagerFactory& factory,
        const AssemblyOptions& options) const
{
    typedef DiscreteVectorValuedLinearOperator<ValueType> DiscreteLinOp;
    typedef DiscreteVectorValuedLinearOperatorSuperposition<ValueType>
            DiscreteSuperposition;

    boost::ptr_vector<DiscreteLinOp> discreteOps;
    for (int i = 0; i < m_terms.size(); ++i)
    {
        std::auto_ptr<DiscreteLinOp> discreteOp =
                m_terms[i].assembleOperator(testPoints, trialSpace, factory, options);
        discreteOps.push_back(discreteOp);
    }
    return std::auto_ptr<DiscreteLinOp>(
                new DiscreteSuperposition(discreteOps, m_multipliers));
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
LinearOperatorSuperposition<ValueType>::assembleWeakForm(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        const IntegrationManagerFactory& factory,
        const AssemblyOptions& options) const
{
    typedef DiscreteScalarValuedLinearOperator<ValueType> DiscreteLinOp;
    typedef DiscreteScalarValuedLinearOperatorSuperposition<ValueType>
            DiscreteSuperposition;

    boost::ptr_vector<DiscreteLinOp> discreteOps;
    for (int i = 0; i < m_terms.size(); ++i)
    {
        std::auto_ptr<DiscreteLinOp> discreteOp =
                m_terms[i].assembleWeakForm(testSpace, trialSpace, factory, options);
        discreteOps.push_back(discreteOp);
    }
    return std::auto_ptr<DiscreteLinOp>(
                new DiscreteSuperposition(discreteOps, m_multipliers));
}

template <typename ValueType>
void LinearOperatorSuperposition<ValueType>::init(
        boost::ptr_vector<LinearOperator<ValueType> >& terms,
        const std::vector<ValueType>& multipliers)
{
    if (terms.size() != multipliers.size())
        throw std::invalid_argument("LinearOperatorSuperposition::init(): "
                                    "incompatible argument lengths");
    if (!terms.empty())
    {
        int testComponentCount_ = terms[0].testComponentCount();
        int trialComponentCount_ = terms[0].trialComponentCount();
        for (int i = 0; i < terms.size(); ++i)
            if (testComponentCount_ != terms[0].testComponentCount() ||
                    trialComponentCount_ != terms[0].trialComponentCount())
                throw std::invalid_argument("LinearOperatorSuperposition::init(): "
                                            "incompatible operator dimensions");
        const int origTermCount = terms.size();
        m_terms.transfer(m_terms.end(), terms);
        const int newTermCount = m_terms.size();
        assert(origTermCount == newTermCount); // documentation of transfer() is
                                               // somewhat difficult to follow...

        m_multipliers = multipliers;
    }
}

#ifdef COMPILE_FOR_FLOAT
template class LinearOperatorSuperposition<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class LinearOperatorSuperposition<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class LinearOperatorSuperposition<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class LinearOperatorSuperposition<std::complex<double> >;
#endif

} // namespace Bempp
