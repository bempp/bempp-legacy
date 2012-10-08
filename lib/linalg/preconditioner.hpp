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

#ifndef bempp_preconditioner_factory_hpp
#define	bempp_preconditioner_factory_hpp

#include "../common/common.hpp"

#include "bempp/common/config_trilinos.hpp"

#ifdef WITH_TRILINOS

#include "Teuchos_RCP.hpp"
#include "Thyra_PreconditionerBase.hpp"
#include "../assembly/discrete_boundary_operator.hpp"
#include "../assembly/discrete_blocked_boundary_operator.hpp"
#include "../assembly/boundary_operator.hpp"
#include "../common/shared_ptr.hpp"
#include "../fiber/scalar_traits.hpp"

#include <vector>

namespace Bempp
{

/** \ingroup linalg
 *  \brief A simple container class to hold pointers to preconditioners.
 *
 */
template<typename ValueType>
class Preconditioner
{
public:
    typedef Teuchos::RCP<const Thyra::PreconditionerBase<ValueType> > TeuchosPreconditionerPtr;
    typedef shared_ptr<const DiscreteBoundaryOperator<ValueType > > DiscreteBoundaryOperatorPtr;
    typedef typename Fiber::ScalarTraits<ValueType>::RealType MagnitudeType;

    Preconditioner(TeuchosPreconditionerPtr precPtr);

    virtual ~Preconditioner();

    /* \brief Return pointer to the actual preconditoner */
    inline const TeuchosPreconditionerPtr& get() const { return m_precPtr;}

private:
    TeuchosPreconditionerPtr m_precPtr;
};

/** \brief Create a preconditioner from a discrete operator.
  */
template<typename ValueType>
Preconditioner<ValueType>
discreteOperatorToPreconditioner(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& discreteOperator);

//template<typename ValueType>
//Preconditioner<ValueType>
//sparseOperatorToPreconditioner(
//        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& discreteOperator);
//
///** \brief Create a preconditioner from a DiscreteBoundaryOperator that
// *  represents a H-Matrix.
// *
// *  Let \f$A\f$ be an H-Matrix. This method computes the H-Matrix LU
// *  factorization \f$A\approx LU\f$ and returns a DiscreteBoundaryOperator
// *  object, which for a given vector \f$y\f$ performs the operation
// *  \f$U^{-1}L^{-1}y\f$.
// *
// *  \param[in] discreteOperator
// *    A DiscreteBoundaryOperator object, which represents an ACA discretized
// *    boundary operator. Only operators that can be * cast to
// *    DiscreteAcaBoundaryOperator are allowed.
// *
// *  \param[in] delta
// *    The accuracy of the approximate H-Matrix LU factorization. A smaller
// *    value means better accuracy, but takes longer to compute (default:
// *    delta=0.01).
// */
//template<typename ValueType>
//Preconditioner<ValueType>
//acaDiscreteOperatorToPreconditioner(
//        const DiscreteBoundaryOperator<ValueType>& discreteOperator,
//        typename Preconditioner<ValueType>::MagnitudeType delta=1E-2);
//
///** \brief Create a block diagonal preconditioner whose blocks are H-Matrix LU
// *  decompositions.
// *
// *  Given operators \f$A_1,\dots,A_N\f$, this method creates a block diagonal
// *  preconditioning operator whose block diagonal entries are solves with the
// *  H-Matrix LU decompositions of \f$A_1,\dots,A_N\f$.
// *
// *  \param[in] opVector
// *
// *     A vector of pointers to DiscreteBoundaryOperator objects, all of which
// *     must be castable to DiscreteAcaBoundaryOperator.
// *
// *  \param[in] deltas
// *     A vector with tolerances for the H-Matrix LU decomposition of each
// *     operator.
// */
//template<typename ValueType>
//Preconditioner<ValueType>
//acaBlockDiagonalPreconditioner(
//        const std::vector<typename Preconditioner<ValueType>::DiscreteBoundaryOperatorPtr>& opVector,
//        const std::vector<typename Preconditioner<ValueType>::MagnitudeType>& deltas);

template<typename ValueType>
Preconditioner<ValueType>
discreteBlockDiagonalPreconditioner(
        const std::vector<shared_ptr<const DiscreteBoundaryOperator<ValueType> > >& opVector);


} // namespace Bempp

#endif /* WITH_TRILINOS */

#endif
