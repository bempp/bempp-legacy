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

namespace Bempp
{

template <typename ValueType>
LinearOperator<ValueType>::LinearOperator(const Space<ValueType>& testSpace,
                                          const Space<ValueType>& trialSpace) :
    m_testSpace(testSpace), m_trialSpace(trialSpace)
{
}

template<typename ValueType>
LinearOperator<ValueType>::LinearOperator(const LinearOperator<ValueType>& other) :
    m_testSpace(other.m_testSpace), m_trialSpace(other.m_trialSpace),
    m_localOperators(other.m_localOperators), m_multipliers(other.m_multipliers)
{
}

template <typename ValueType>
LinearOperator<ValueType>::~LinearOperator()
{
}

template <typename ValueType>
void LinearOperator<ValueType>::assemble(const LocalAssemblerFactory& factory,
                                         const AssemblyOptions& options)
{
    m_discreteOperator = this->assembleWeakForm(factory, options);
}

template <typename ValueType>
bool LinearOperator<ValueType>::isAssembled() const
{
    return m_discreteOperator.get() != NULL;
}

template <typename ValueType>
void LinearOperator<ValueType>::apply(
        const TranspositionMode trans,
        const GridFunction<ValueType>& x_in,
        GridFunction<ValueType>& y_inout,
        ValueType alpha, ValueType beta) const
{
    if (!this->isAssembled())
        throw std::runtime_error("LinearOperator::apply(): "
                                 "operator is not assembled");

    // Sanity test
    if (&m_trialSpace != &x_in.space() || &m_testSpace != &y_inout.space())
        throw std::runtime_error("LinearOperator::apply(): Spaces don't match");

    // Extract coefficient vectors
    arma::Col<ValueType> xVals = x_in.coefficients().asArmadilloVector();
    arma::Col<ValueType> yVals = y_inout.coefficients().asArmadilloVector();

    // Apply operator and assign the result to y_inout's coefficients
    m_discreteOperator->apply(trans, xVals, yVals, alpha, beta);
    // Note: all these armadillo<->vector conversions are horribly
    // inefficient and unnecessary. But in order to get rid of them, we need
    // first to make interfaces to the Trilinos and fallback
    // DiscreteLinearOperator::apply() compatible.
    // Perhaps by declaring an asPtrToBaseVector method in Vector...
    y_inout.setCoefficients(Vector<ValueType>(yVals));
}

template <typename ValueType>
const DiscreteLinearOperator<ValueType>&
LinearOperator<ValueType>::assembledDiscreteLinearOperator() const
{
    if (!isAssembled())
        throw std::runtime_error("LinearOperator::assembledDiscreteLinearOperator(): "
                                 "operator is not assembled");
    return *m_discreteOperator;
}

template <typename ValueType>
const std::vector<ElementaryLinearOperator<ValueType> const*>&
LinearOperator<ValueType>::localOperators() const
{
    return m_localOperators;
}

template <typename ValueType>
const std::vector<ValueType >& LinearOperator<ValueType>::multipliers() const
{
    return m_multipliers;
}

template <typename ValueType>
const Space<ValueType>& LinearOperator<ValueType>::testSpace() const
{
    return m_testSpace;
}

template <typename ValueType>
const Space<ValueType>& LinearOperator<ValueType>::trialSpace() const
{
    return m_trialSpace;
}

template <typename ValueType>
void LinearOperator<ValueType>::addLocalOperatorsAndMultipliers(
        const std::vector<ElementaryLinearOperator<ValueType> const*>& localOperators,
        const std::vector<ValueType>& multipliers)
{
    m_localOperators.insert(m_localOperators.end(),
                            localOperators.begin(), localOperators.end());
    m_multipliers.insert(m_multipliers.end(),
                         multipliers.begin(), multipliers.end());
}

template <typename ValueType>
LinearOperatorSuperposition<ValueType> operator+(
        const LinearOperator<ValueType>& op1,
        const LinearOperator<ValueType>& op2)
{
    return LinearOperatorSuperposition<ValueType>(op1, op2);
}

template <typename ValueType>
LinearOperatorSuperposition<ValueType> operator-(
        const LinearOperator<ValueType>& op1,
        const LinearOperator<ValueType>& op2)
{
    return op1 + (-1. * op2);
}

template <typename ValueType>
LinearOperatorSuperposition<ValueType> operator*(
        const LinearOperator<ValueType>& op, const ValueType& scalar)
{
    return LinearOperatorSuperposition<ValueType>(op, scalar);
}

template <typename ValueType>
LinearOperatorSuperposition<ValueType> operator*(
        const ValueType& scalar, const LinearOperator<ValueType>& op)
{
    return LinearOperatorSuperposition<ValueType>(op, scalar);
}

template <typename ValueType>
LinearOperatorSuperposition<ValueType> operator/(
        const LinearOperator<ValueType>& op, const ValueType& scalar)
{
    if (scalar == 0)
        throw std::runtime_error("LinearOperatorSuperposition::operator/(): "
                                 "Division by zero");

    return LinearOperatorSuperposition<ValueType>(op, 1. / scalar);
}

template <typename ValueType>
GridFunction<ValueType> operator*(
        const LinearOperator<ValueType>& op,
        const GridFunction<ValueType>& fun)
{
    const Space<ValueType>& space = op.testSpace();
    arma::Col<ValueType> coeffs = arma::zeros(space.globalDofCount());
    GridFunction<ValueType> result(space, coeffs);
    op.apply(NO_TRANSPOSE, fun, result, 1., 0.);
    return result;
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
