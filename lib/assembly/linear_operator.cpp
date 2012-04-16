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

#include "linear_operator.hpp"
#include "linear_operator_superposition.hpp"
#include "discrete_linear_operator.hpp"
#include <stdexcept>

namespace Bempp{

template <typename ValueType>
LinearOperator<ValueType>::LinearOperator(const Space<ValueType>& testSpace, const Space<ValueType>& trialSpace) :
    m_testSpace(testSpace), m_trialSpace(trialSpace) {}

template <typename ValueType>
void LinearOperator<ValueType>::assemble(const LocalAssemblerFactory& factory, const AssemblyOptions& options){
    m_discreteOperator=this->assembleWeakForm(factory,options);
}

template <typename ValueType>
bool LinearOperator<ValueType>::isAssembled() const{
    return m_discreteOperator.get()!=NULL;
}

template <typename ValueType>
GridFunction<ValueType> LinearOperator<ValueType>::apply(
        const GridFunction<ValueType>& x, const GridFunction<ValueType>&  y, ValueType alpha, ValueType beta) const{
    if (!this->isAssembled()) throw std::runtime_error("Operator is not assembled");

    // Sanity test
    if ((m_trialSpace!=x.space()) || (m_testSpace!=y.space())) throw std::runtime_error("Spaces don't match");

    // Extract coefficient vectors

    arma::Col<ValueType> xVals=x.coefficients().asArmadilloVector();
    arma::Col<ValueType> yVals=y.coefficients().asArmadilloVector();

    // Apply operator and return new GridFunction

    m_discreteOperator->apply(Bempp::NO_TRANSPOSE,xVals,yVals,alpha,beta);
    return GridFunction<ValueType>(m_testSpace,yVals);

}

template <typename ValueType>
GridFunction<ValueType> LinearOperator<ValueType>::apply(const GridFunction<ValueType>& x) const{
     if (!this->isAssembled()) throw std::runtime_error("Operator is not assembled");

     //Sanity test

     if (m_trialSpace!=x.space()) throw std::runtime_error("Spaces don't match");

     arma::Col<ValueType> xVals=x.coefficients().asArmadilloVector();
     arma::Col<ValueType> yVals=arma::zeros(m_testSpace.globalDofCount());
     m_discreteOperator->apply(Bempp::NO_TRANSPOSE,xVals,yVals,1.0,1.0);
     return GridFunction<ValueType>(m_testSpace,yVals);

}

template <typename ValueType>
const DiscreteLinearOperator<ValueType>& LinearOperator<ValueType>::getDiscreteLinearOperator() const{
    return *m_discreteOperator;
}

template<typename ValueType>
LinearOperator<ValueType>::LinearOperator(const LinearOperator<ValueType>& linOp) :
    m_testSpace(linOp.m_testSpace), m_trialSpace(linOp.m_trialSpace),
    m_localOperators(linOp.m_localOperators), m_multipliers(linOp.m_multipliers){}

template <typename ValueType>
LinearOperator<ValueType>::~LinearOperator(){}

template <typename ValueType>
const std::vector<ElementaryLinearOperator<ValueType> const* >& LinearOperator<ValueType>::getLocalOperators() const {
    return m_localOperators;
}

template <typename ValueType>
const std::vector<ValueType >& LinearOperator<ValueType>::getMultipliers() const {
    return m_multipliers;
}

template <typename ValueType>
const Space<ValueType>& LinearOperator<ValueType>::getTestSpace() const{
    return m_testSpace;
}

template <typename ValueType>
const Space<ValueType>& LinearOperator<ValueType>::getTrialSpace() const{
    return m_trialSpace;
}

template <typename ValueType>
bool LinearOperator<ValueType>::checkSpaces(const LinearOperator<ValueType>& linOp) const{
    return ((m_testSpace==linOp.getTestSpace()) && (m_trialSpace==linOp.getTrialSpace()) );
}


template <typename ValueType>
void LinearOperator<ValueType>::addLocalOperatorsMultipliers(const std::vector<ElementaryLinearOperator<ValueType> const*>& localOperators,
                                                             const std::vector<ValueType>& multipliers)
{
    m_localOperators.insert(m_localOperators.end(),localOperators.begin(),localOperators.end());
    m_multipliers.insert(m_multipliers.end(),multipliers.begin(),multipliers.end());

}

template <typename ValueType>
LinearOperatorSuperposition<ValueType> operator+(const LinearOperator<ValueType>& op1, const LinearOperator<ValueType>& op2){
    return LinearOperatorSuperposition<ValueType>(op1,op2);
}

template <typename ValueType>
LinearOperatorSuperposition<ValueType> operator-(const LinearOperator<ValueType>& op1, const LinearOperator<ValueType>& op2){
    return op1+(-1.*op2);
}


template <typename ValueType>
LinearOperatorSuperposition<ValueType> operator*(const LinearOperator<ValueType>& op, const ValueType& scalar){
    return LinearOperatorSuperposition<ValueType>(op,scalar);
}

template <typename ValueType>
LinearOperatorSuperposition<ValueType> operator*(const ValueType& scalar, const LinearOperator<ValueType>& op){
    return LinearOperatorSuperposition<ValueType>(op,scalar);
}

template <typename ValueType>
LinearOperatorSuperposition<ValueType> operator/(const LinearOperator<ValueType>& op, const ValueType& scalar){

    if (scalar==0) throw std::runtime_error("Division by zero");

    return LinearOperatorSuperposition<ValueType>(op,1./scalar);
}

template <typename ValueType>
GridFunction<ValueType> operator*(const LinearOperator<ValueType>& op, const GridFunction<ValueType>& fun){
    return op.apply(fun);
}


#ifdef COMPILE_FOR_FLOAT
template class LinearOperator<float>;
template LinearOperatorSuperposition<float> operator+(const LinearOperator<float>& op1, const LinearOperator<float>& op2);
template LinearOperatorSuperposition<float> operator-(const LinearOperator<float>& op1, const LinearOperator<float>& op2);
template LinearOperatorSuperposition<float> operator*(const LinearOperator<float>& op, const float& scalar);
template LinearOperatorSuperposition<float> operator*(const float& scalar, const LinearOperator<float>& op);
template LinearOperatorSuperposition<float> operator/(const LinearOperator<float>& op, const float& scalar);
template GridFunction<float> operator*(const LinearOperator<float>& op, const GridFunction<float>& fun);
#endif
#ifdef COMPILE_FOR_DOUBLE
template class LinearOperator<double>;
template LinearOperatorSuperposition<double> operator+(const LinearOperator<double>& op1, const LinearOperator<double>& op2);
template LinearOperatorSuperposition<double> operator-(const LinearOperator<double>& op1, const LinearOperator<double>& op2);
template LinearOperatorSuperposition<double> operator*(const LinearOperator<double>& op, const double& scalar);
template LinearOperatorSuperposition<double> operator*(const double& scalar, const LinearOperator<double>& op);
template LinearOperatorSuperposition<double> operator/(const LinearOperator<double>& op, const double& scalar);
template GridFunction<double> operator*(const LinearOperator<double>& op, const GridFunction<double>& fun);
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class LinearOperator<std::complex<float> >;
template LinearOperatorSuperposition<std::complex<float> > operator+(const LinearOperator<std::complex<float> >& op1, const LinearOperator<std::complex<float> >& op2);
template LinearOperatorSuperposition<std::complex<float> > operator-(const LinearOperator<std::complex<float> >& op1, const LinearOperator<std::complex<float> >& op2);
template LinearOperatorSuperposition<std::complex<float> > operator*(const LinearOperator<std::complex<float> >& op, const std::complex<float>& scalar);
template LinearOperatorSuperposition<std::complex<float> > operator*(const std::complex<float>& scalar, const LinearOperator<fstd::complex<float> >& op);
template LinearOperatorSuperposition<std::complex<float> > operator/(const LinearOperator<std::complex<float> >& op, const std::complex<float>& scalar);
template GridFunction<std::complex<float> > operator*(const LinearOperator<std::complex<float> >& op, const GridFunction<std::complex<float> >& fun);
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class LinearOperator<std::complex<double> >;
template LinearOperatorSuperposition<std::complex<double> > operator+(const LinearOperator<std::complex<double> >& op1, const LinearOperator<std::complex<double> >& op2);
template LinearOperatorSuperposition<std::complex<double> > operator-(const LinearOperator<std::complex<double> >& op1, const LinearOperator<std::complex<double> >& op2);
template LinearOperatorSuperposition<std::complex<double> > operator*(const LinearOperator<std::complex<double> >& op, const std::complex<double>& scalar);
template LinearOperatorSuperposition<std::complex<double> > operator*(const std::complex<double>& scalar, const LinearOperator<std::complex<double> >& op);
template LinearOperatorSuperposition<std::complex<double> > operator/(const LinearOperator<std::complex<double> >& op, const std::complex<double>& scalar);
template GridFunction<std::complex<double> > operator*(const LinearOperator<std::complex<double> >& op, const GridFunction<std::complex<double> >& fun);
#endif

}
