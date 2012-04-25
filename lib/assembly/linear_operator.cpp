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
#include "../fiber/explicit_instantiation.hpp"

#include <stdexcept>

namespace Bempp
{

template <typename ArgumentType, typename ResultType>
LinearOperator<ArgumentType, ResultType>::
LinearOperator(const Space<ArgumentType>& testSpace,
               const Space<ArgumentType>& trialSpace) :
    m_testSpace(testSpace), m_trialSpace(trialSpace)
{
}

template <typename ArgumentType, typename ResultType>
LinearOperator<ArgumentType, ResultType>::LinearOperator(
        const LinearOperator<ArgumentType, ResultType>& other) :
    m_testSpace(other.m_testSpace), m_trialSpace(other.m_trialSpace),
    m_localOperators(other.m_localOperators), m_multipliers(other.m_multipliers)
{
}

template <typename ArgumentType, typename ResultType>
LinearOperator<ArgumentType, ResultType>::~LinearOperator()
{
}

template <typename ArgumentType, typename ResultType>
void LinearOperator<ArgumentType, ResultType>::assemble(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options)
{
    m_discreteOperator = this->assembleWeakForm(factory, options);
}

template <typename ArgumentType, typename ResultType>
bool LinearOperator<ArgumentType, ResultType>::isAssembled() const
{
    return m_discreteOperator.get() != NULL;
}

template <typename ArgumentType, typename ResultType>
void LinearOperator<ArgumentType, ResultType>::apply(
        const TranspositionMode trans,
        const GridFunction<ArgumentType, ResultType>& x_in,
        GridFunction<ArgumentType, ResultType>& y_inout,
        ResultType alpha, ResultType beta) const
{
    if (!this->isAssembled())
        throw std::runtime_error("LinearOperator::apply(): "
                                 "operator is not assembled");

    // Sanity test
    if (&m_trialSpace != &x_in.space() || &m_testSpace != &y_inout.space())
        throw std::runtime_error("LinearOperator::apply(): Spaces don't match");

    // Extract coefficient vectors
    arma::Col<ResultType> xVals = x_in.coefficients().asArmadilloVector();
    arma::Col<ResultType> yVals = y_inout.coefficients().asArmadilloVector();

    // Apply operator and assign the result to y_inout's coefficients
    m_discreteOperator->apply(trans, xVals, yVals, alpha, beta);
    // Note: all these armadillo<->vector conversions are horribly
    // inefficient and unnecessary. But in order to get rid of them, we need
    // first to make interfaces to the Trilinos and fallback
    // DiscreteLinearOperator::apply() compatible.
    // Perhaps by declaring an asPtrToBaseVector method in Vector...
    y_inout.setCoefficients(Vector<ResultType>(yVals));
}

template <typename ArgumentType, typename ResultType>
const DiscreteLinearOperator<ResultType>&
LinearOperator<ArgumentType, ResultType>::assembledDiscreteLinearOperator() const
{
    if (!isAssembled())
        throw std::runtime_error("LinearOperator::assembledDiscreteLinearOperator(): "
                                 "operator is not assembled");
    return *m_discreteOperator;
}

template <typename ArgumentType, typename ResultType>
const std::vector<ElementaryLinearOperator<ArgumentType, ResultType> const*>&
LinearOperator<ArgumentType, ResultType>::localOperators() const
{
    return m_localOperators;
}

template <typename ArgumentType, typename ResultType>
const std::vector<ResultType>& LinearOperator<ArgumentType, ResultType>::multipliers() const
{
    return m_multipliers;
}

template <typename ArgumentType, typename ResultType>
const Space<ArgumentType>& LinearOperator<ArgumentType, ResultType>::testSpace() const
{
    return m_testSpace;
}

template <typename ArgumentType, typename ResultType>
const Space<ArgumentType>& LinearOperator<ArgumentType, ResultType>::trialSpace() const
{
    return m_trialSpace;
}

template <typename ArgumentType, typename ResultType>
void LinearOperator<ArgumentType, ResultType>::addLocalOperatorsAndMultipliers(
        const std::vector<ElementaryLinearOperator<ArgumentType, ResultType> const*>&
        localOperators,
        const std::vector<ResultType>& multipliers)
{
    m_localOperators.insert(m_localOperators.end(),
                            localOperators.begin(), localOperators.end());
    m_multipliers.insert(m_multipliers.end(),
                         multipliers.begin(), multipliers.end());
}

template <typename ArgumentType, typename ResultType>
LinearOperatorSuperposition<ArgumentType, ResultType> operator+(
        const LinearOperator<ArgumentType, ResultType>& op1,
        const LinearOperator<ArgumentType, ResultType>& op2)
{
    return LinearOperatorSuperposition<ArgumentType, ResultType>(op1, op2);
}

template <typename ArgumentType, typename ResultType>
LinearOperatorSuperposition<ArgumentType, ResultType> operator-(
        const LinearOperator<ArgumentType, ResultType>& op1,
        const LinearOperator<ArgumentType, ResultType>& op2)
{
    return op1 + (static_cast<ResultType>(-1.) * op2);
}

template <typename ArgumentType, typename ResultType, typename ScalarType>
LinearOperatorSuperposition<ArgumentType, ResultType> operator*(
        const LinearOperator<ArgumentType, ResultType>& op, const ScalarType& scalar)
{
    return LinearOperatorSuperposition<ArgumentType, ResultType>(op, scalar);
}

template <typename ArgumentType, typename ResultType, typename ScalarType>
LinearOperatorSuperposition<ArgumentType, ResultType> operator*(
        const ScalarType& scalar, const LinearOperator<ArgumentType, ResultType>& op)
{
    return LinearOperatorSuperposition<ArgumentType, ResultType>(op, scalar);
}

template <typename ArgumentType, typename ResultType, typename ScalarType>
LinearOperatorSuperposition<ArgumentType, ResultType> operator/(
        const LinearOperator<ArgumentType, ResultType>& op, const ScalarType& scalar)
{
    if (scalar == 0.)
        throw std::runtime_error("LinearOperatorSuperposition::operator/(): "
                                 "Division by zero");

    return LinearOperatorSuperposition<ArgumentType, ResultType>(
                op, static_cast<ResultType>(1.) / scalar);
}

template <typename ArgumentType, typename ResultType>
GridFunction<ArgumentType, ResultType> operator*(
        const LinearOperator<ArgumentType, ResultType>& op,
        const GridFunction<ArgumentType, ResultType>& fun)
{
    const Space<ArgumentType>& space = op.testSpace();
    arma::Col<ResultType> coeffs(space.globalDofCount());
    coeffs.fill(0.);
    GridFunction<ArgumentType, ResultType> result(space, coeffs);
    op.apply(NO_TRANSPOSE, fun, result, 1., 0.);
    return result;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(LinearOperator);

//#ifdef COMPILE_FOR_FLOAT
//template class LinearOperator<float>;
//template LinearOperatorSuperposition<float> operator+(
//    const LinearOperator<float>& op1, const LinearOperator<float>& op2);
//template LinearOperatorSuperposition<float> operator-(
//    const LinearOperator<float>& op1, const LinearOperator<float>& op2);
//template LinearOperatorSuperposition<float> operator*(
//    const LinearOperator<float>& op, const float& scalar);
//template LinearOperatorSuperposition<float> operator*(
//    const LinearOperator<float>& op, const double& scalar);
//template LinearOperatorSuperposition<float> operator*(
//    const float& scalar, const LinearOperator<float>& op);
//template LinearOperatorSuperposition<float> operator*(
//    const double& scalar, const LinearOperator<float>& op);
//template LinearOperatorSuperposition<float> operator/(
//    const LinearOperator<float>& op, const float& scalar);
//template LinearOperatorSuperposition<float> operator/(
//    const LinearOperator<float>& op, const double& scalar);
//template GridFunction<float> operator*(
//    const LinearOperator<float>& op, const GridFunction<float>& fun);
//#endif
//#ifdef COMPILE_FOR_DOUBLE
//template class LinearOperator<double>;
//template LinearOperatorSuperposition<double> operator+(
//    const LinearOperator<double>& op1, const LinearOperator<double>& op2);
//template LinearOperatorSuperposition<double> operator-(
//    const LinearOperator<double>& op1, const LinearOperator<double>& op2);
//template LinearOperatorSuperposition<double> operator*(
//    const LinearOperator<double>& op, const float& scalar);
//template LinearOperatorSuperposition<double> operator*(
//    const LinearOperator<double>& op, const double& scalar);
//template LinearOperatorSuperposition<double> operator*(
//    const float& scalar, const LinearOperator<double>& op);
//template LinearOperatorSuperposition<double> operator*(
//    const double& scalar, const LinearOperator<double>& op);
//template LinearOperatorSuperposition<double> operator/(
//    const LinearOperator<double>& op, const float& scalar);
//template LinearOperatorSuperposition<double> operator/(
//    const LinearOperator<double>& op, const double& scalar);
//template GridFunction<double> operator*(
//    const LinearOperator<double>& op, const GridFunction<double>& fun);
//#endif
//#ifdef COMPILE_FOR_COMPLEX_FLOAT
//#include <complex>
//template class LinearOperator<std::complex<float> >;
//template LinearOperatorSuperposition<std::complex<float> > operator+(
//    const LinearOperator<std::complex<float> >& op1,
//    const LinearOperator<std::complex<float> >& op2);
//template LinearOperatorSuperposition<std::complex<float> > operator-(
//    const LinearOperator<std::complex<float> >& op1,
//    const LinearOperator<std::complex<float> >& op2);

//template LinearOperatorSuperposition<std::complex<float> > operator*(
//    const LinearOperator<std::complex<float> >& op,
//    const float& scalar);
//template LinearOperatorSuperposition<std::complex<float> > operator*(
//    const LinearOperator<std::complex<float> >& op,
//    const double& scalar);
//template LinearOperatorSuperposition<std::complex<float> > operator*(
//    const LinearOperator<std::complex<float> >& op,
//    const std::complex<float>& scalar);
//template LinearOperatorSuperposition<std::complex<float> > operator*(
//    const LinearOperator<std::complex<float> >& op,
//    const std::complex<double>& scalar);

//template LinearOperatorSuperposition<std::complex<float> > operator*(
//    const float& scalar,
//    const LinearOperator<fstd::complex<float> >& op);
//template LinearOperatorSuperposition<std::complex<float> > operator*(
//    const double& scalar,
//    const LinearOperator<fstd::complex<float> >& op);
//template LinearOperatorSuperposition<std::complex<float> > operator*(
//    const std::complex<float>& scalar,
//    const LinearOperator<fstd::complex<float> >& op);
//template LinearOperatorSuperposition<std::complex<float> > operator*(
//    const std::complex<double>& scalar,
//    const LinearOperator<fstd::complex<float> >& op);

//template LinearOperatorSuperposition<std::complex<float> > operator/(
//    const LinearOperator<std::complex<float> >& op,
//    const float& scalar);
//template LinearOperatorSuperposition<std::complex<float> > operator/(
//    const LinearOperator<std::complex<float> >& op,
//    const double& scalar);
//template LinearOperatorSuperposition<std::complex<float> > operator/(
//    const LinearOperator<std::complex<float> >& op,
//    const std::complex<float>& scalar);
//template LinearOperatorSuperposition<std::complex<float> > operator/(
//    const LinearOperator<std::complex<float> >& op,
//    const std::complex<double>& scalar);

//template GridFunction<std::complex<float> > operator*(
//    const LinearOperator<std::complex<float> >& op,
//    const GridFunction<std::complex<float> >& fun);
//bc;
//#endif
//#ifdef COMPILE_FOR_COMPLEX_DOUBLE
//#include <complex>
//template class LinearOperator<std::complex<double> >;
//template LinearOperatorSuperposition<std::complex<double> > operator+(
//    const LinearOperator<std::complex<double> >& op1,
//    const LinearOperator<std::complex<double> >& op2);
//template LinearOperatorSuperposition<std::complex<double> > operator-(
//    const LinearOperator<std::complex<double> >& op1,
//    const LinearOperator<std::complex<double> >& op2);

//template LinearOperatorSuperposition<std::complex<double> > operator*(
//    const LinearOperator<std::complex<double> >& op,
//    const float& scalar);
//template LinearOperatorSuperposition<std::complex<double> > operator*(
//    const LinearOperator<std::complex<double> >& op,
//    const double& scalar);
//template LinearOperatorSuperposition<std::complex<double> > operator*(
//    const LinearOperator<std::complex<double> >& op,
//    const std::complex<float>& scalar);
//template LinearOperatorSuperposition<std::complex<double> > operator*(
//    const LinearOperator<std::complex<double> >& op,
//    const std::complex<double>& scalar);

//template LinearOperatorSuperposition<std::complex<double> > operator*(
//    const float& scalar,
//    const LinearOperator<fstd::complex<double> >& op);
//template LinearOperatorSuperposition<std::complex<double> > operator*(
//    const double& scalar,
//    const LinearOperator<fstd::complex<double> >& op);
//template LinearOperatorSuperposition<std::complex<double> > operator*(
//    const std::complex<float>& scalar,
//    const LinearOperator<fstd::complex<double> >& op);
//template LinearOperatorSuperposition<std::complex<double> > operator*(
//    const std::complex<double>& scalar,
//    const LinearOperator<fstd::complex<double> >& op);

//template LinearOperatorSuperposition<std::complex<double> > operator/(
//    const LinearOperator<std::complex<double> >& op,
//    const float& scalar);
//template LinearOperatorSuperposition<std::complex<double> > operator/(
//    const LinearOperator<std::complex<double> >& op,
//    const double& scalar);
//template LinearOperatorSuperposition<std::complex<double> > operator/(
//    const LinearOperator<std::complex<double> >& op,
//    const std::complex<float>& scalar);
//template LinearOperatorSuperposition<std::complex<double> > operator/(
//    const LinearOperator<std::complex<double> >& op,
//    const std::complex<double>& scalar);

//template GridFunction<std::complex<double> > operator*(
//    const LinearOperator<std::complex<double> >& op,
//    const GridFunction<std::complex<double> >& fun);
//#endif

}
