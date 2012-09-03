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

#include "bempp/common/config_trilinos.hpp"
#include "bempp/common/config_ahmed.hpp"

#if defined(WITH_TRILINOS) && defined(WITH_AHMED)

#include "preconditioner.hpp"
#include "../assembly/discrete_aca_boundary_operator.hpp"
#include "../assembly/discrete_boundary_operator.hpp"
#include "../assembly/aca_approximate_lu_inverse.hpp"
#include "../assembly/discrete_aca_boundary_operator.hpp"
#include "../assembly/aca_approximate_lu_inverse.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/_2d_array.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../assembly/discrete_blocked_boundary_operator.hpp"

#include <Teuchos_RCP.hpp>
#include <Thyra_LinearOpBase.hpp>
#include <Thyra_PreconditionerBase.hpp>
#include <Thyra_DefaultPreconditioner.hpp>

namespace Bempp
{

template<typename ValueType>
Preconditioner<ValueType>::Preconditioner(TeuchosPreconditionerPtr precPtr):
    m_precPtr(precPtr){
}

template<typename ValueType>
Preconditioner<ValueType>::~Preconditioner(){}


template<typename ValueType>
Preconditioner<ValueType>
acaDiscreteOperatorToPreconditioner(
        const DiscreteBoundaryOperator<ValueType>& discreteOperator,
        typename Preconditioner<ValueType>::MagnitudeType delta)
{
    const DiscreteAcaBoundaryOperator<ValueType>& discreteAcaOperator =
                DiscreteAcaBoundaryOperator<ValueType>::castToAca(discreteOperator);

    Teuchos::RCP<const Thyra::LinearOpBase<ValueType> > op(
                new AcaApproximateLuInverse<ValueType>(discreteAcaOperator,delta));

    typename Preconditioner<ValueType>::TeuchosPreconditionerPtr precOp =
            Teuchos::rcp_static_cast<const Thyra::PreconditionerBase<ValueType> >(
                Thyra::unspecifiedPrec(op));

    return Preconditioner<ValueType>(precOp);
}

template<typename ValueType>
Preconditioner<ValueType>
acaBlockDiagonalPreconditioner(
        const std::vector<typename Preconditioner<ValueType>::DiscreteBoundaryOperatorPtr>& opVector,
        const std::vector<typename Preconditioner<ValueType>::MagnitudeType>& deltas)
{
    size_t n = opVector.size();
    if (n != deltas.size())
        throw std::runtime_error("acaBlockDiagonalPreconditioner: "
                                 "Input arrays must have same length");
    if (n == 0)
        throw std::runtime_error("acaBlockDiagonalPreconditioner: "
                                 "Input arrays must not be empty");

    Fiber::_2dArray<typename Preconditioner<ValueType>::DiscreteBoundaryOperatorPtr> opStructure(n,n);
    std::vector<size_t> rowCounts;
    std::vector<size_t> columnCounts;
    for (size_t i=0;i++;i<n){
        opStructure(i,i)= typename Preconditioner<ValueType>::DiscreteBoundaryOperatorPtr(
                    new AcaApproximateLuInverse<ValueType>(
                        DiscreteAcaBoundaryOperator<ValueType>::castToAca(*(opVector[i]))
                        ,deltas[i]));
        rowCounts[i]=opVector[i]->rowCount();
        columnCounts[i]=opVector[i]->columnCount();
    }

    Teuchos::RCP<const Thyra::LinearOpBase<ValueType> > op
            (new DiscreteBlockedBoundaryOperator<ValueType>(
                                            opStructure,
                                            rowCounts,
                                            columnCounts));

    typename Preconditioner<ValueType>::TeuchosPreconditionerPtr precOp =
            Teuchos::rcp_static_cast<const Thyra::PreconditionerBase<ValueType> >(
                Thyra::unspecifiedPrec(op));
    return Preconditioner<ValueType>(precOp);

}

#define INSTANTIATE_FREE_FUNCTIONS( VALUE ) \
    template Preconditioner< VALUE > \
    acaDiscreteOperatorToPreconditioner(const DiscreteBoundaryOperator< VALUE >& discreteOperator, \
                                        Preconditioner< VALUE >::MagnitudeType delta=1E-2); \
    template Preconditioner< VALUE > \
    acaBlockDiagonalPreconditioner(const std::vector<Preconditioner< VALUE >::DiscreteBoundaryOperatorPtr>& opVector, \
                                   const std::vector<Preconditioner< VALUE >::MagnitudeType>& deltas);




FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(Preconditioner);
FIBER_ITERATE_OVER_VALUE_TYPES(INSTANTIATE_FREE_FUNCTIONS);


} // namespace Bempp

#endif // WITH_TRILINOS && WITH_AHMED
