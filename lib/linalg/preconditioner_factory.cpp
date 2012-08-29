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

#include "preconditioner_factory.hpp"
#include "../assembly/discrete_aca_boundary_operator.hpp"
#include "../assembly/discrete_boundary_operator.hpp"
#include "../assembly/aca_approximate_lu_inverse.hpp"
#include "../assembly/discrete_aca_boundary_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <Teuchos_RCP.hpp>
#include <Thyra_LinearOpBase.hpp>
#include <Thyra_PreconditionerBase.hpp>
#include <Thyra_DefaultPreconditioner.hpp>

namespace Bempp
{

template<typename ValueType>
Teuchos::RCP<const Thyra::PreconditionerBase<ValueType> >
PreconditionerFactory<ValueType>::acaOperatorToPreconditioner(
        const DiscreteBoundaryOperator<ValueType>& discreteOperator, const double delta)
{
    const DiscreteAcaBoundaryOperator<ValueType>& discreteAcaOperator =
                DiscreteAcaBoundaryOperator<ValueType>::castToAca(discreteOperator);

    Teuchos::RCP<const Thyra::LinearOpBase<ValueType> > precOp(
                new AcaApproximateLuInverse<ValueType>(discreteAcaOperator,delta));

    Teuchos::RCP<const Thyra::PreconditionerBase<ValueType> > preconditioner =
            Teuchos::rcp_static_cast<const Thyra::PreconditionerBase<ValueType> >(
                Thyra::unspecifiedPrec(precOp));

    return preconditioner;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(PreconditionerFactory);

} // namespace Bempp

#endif // WITH_TRILINOS && WITH_AHMED
