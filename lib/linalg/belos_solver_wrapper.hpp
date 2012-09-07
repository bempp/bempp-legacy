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

#ifndef belos_solver_wrapper_hpp
#define belos_solver_wrapper_hpp

#include "../common/common.hpp"

#include "belos_solver_wrapper_fwd.hpp"
#include "../common/scalar_traits.hpp"

namespace Thyra
{
/** \cond FORWARD_DECL */
template <typename ValueType> class PreconditionerBase;
template <typename ValueType> class LinearOpWithSolveBase;
/** \endcond */
}

namespace Bempp
{

template <typename ValueType>
class BelosSolverWrapper
{
public:
    typedef typename ScalarTraits<ValueType>::RealType MagnitudeType;

    BelosSolverWrapper(
            const Teuchos::RCP<const Thyra::LinearOpBase<ValueType> >& linOp);

    ~BelosSolverWrapper();

    void setPreconditioner(
            const Teuchos::RCP<const Thyra::PreconditionerBase<ValueType> >& preconditioner);

    void initializeSolver(
            const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    Thyra::SolveStatus<MagnitudeType> solve(
            const Thyra::EOpTransp trans,
            const Thyra::MultiVectorBase<ValueType>& rhs,
            const Teuchos::Ptr<Thyra::MultiVectorBase<ValueType> >& sol) const;

private:
    Teuchos::RCP<const Thyra::LinearOpBase<ValueType> > m_linOp;
    Teuchos::RCP<const Thyra::PreconditionerBase<ValueType> > m_preconditioner;
    Teuchos::RCP<const Thyra::LinearOpWithSolveBase<MagnitudeType> > m_linOpWithSolve;
};

} // namespace Bempp

#endif
