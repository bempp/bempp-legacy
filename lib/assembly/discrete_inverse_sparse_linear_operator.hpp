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

#ifndef bempp_discrete_inverse_sparse_linear_operator_hpp
#define bempp_discrete_inverse_sparse_linear_operator_hpp

#include "config_trilinos.hpp"

#ifdef WITH_TRILINOS
#include "discrete_linear_operator.hpp"

#include <memory>

#include <Teuchos_RCP.hpp>
#include <Thyra_SpmdVectorSpaceBase_decl.hpp>
class Amesos_BaseSolver;
class Epetra_LinearProblem;
class Epetra_CrsMatrix;

namespace Bempp
{

/** \ingroup assembly
 *  \brief Discrete linear operator being an inverse of a sparse matrix.
 */
template <typename ValueType>
class DiscreteInverseSparseLinearOperator :
        public DiscreteLinearOperator<ValueType>
{
public:
    DiscreteInverseSparseLinearOperator(const Epetra_CrsMatrix& mat,
                                        bool m_symmetric = false);

    virtual void dump() const;

    virtual arma::Mat<ValueType> asMatrix() const;

    virtual unsigned int rowCount() const;
    virtual unsigned int columnCount() const;

    virtual void addBlock(const std::vector<int>& rows,
                          const std::vector<int>& cols,
                          const ValueType alpha,
                          arma::Mat<ValueType>& block) const;

public:
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > domain() const;
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > range() const;

protected:
    virtual bool opSupportedImpl(Thyra::EOpTransp M_trans) const;
    virtual void applyImpl(
            const Thyra::EOpTransp M_trans,
            const Thyra::MultiVectorBase<ValueType> &X_in,
            const Teuchos::Ptr<Thyra::MultiVectorBase<ValueType> > &Y_inout,
            const ValueType alpha,
            const ValueType beta) const;

private:
    virtual void applyBuiltInImpl(const TranspositionMode trans,
                                  const arma::Col<ValueType>& x_in,
                                  arma::Col<ValueType>& y_inout,
                                  const ValueType alpha,
                                  const ValueType beta) const;

private:
    std::auto_ptr<Epetra_LinearProblem> m_problem;
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType> > m_space;
    bool m_symmetric;
    std::auto_ptr<Amesos_BaseSolver> m_solver;
};

} // namespace Bempp

#endif // WITH_TRILINOS

#endif
