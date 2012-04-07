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

#include "../common/config_trilinos.hpp"

#ifndef aca_preconditioner_factory_hpp
#define	aca_preconditioner_factory_hpp

#ifdef WITH_TRILINOS

#include <Teuchos_RCP.hpp>

#include "../assembly/discrete_scalar_valued_linear_operator.hpp"


namespace Bempp {

template<typename ValueType>
class AcaPreconditionerFactory {
    
    public:
        typedef Teuchos::RCP<const Thyra::PreconditionerBase<ValueType> > ThyraPreconditioner;
        
        static ThyraPreconditioner getThyraPreconditioner
                (DiscreteScalarValuedLinearOperator<ValueType>& discreteOperator, const double delta=0.1);
        
        
};


}

#endif /* WITH_TRILINOS */

#endif	/* aca_preconditioner_factory_hpp */

//       RCP<const Thyra::LinearOpBase<double> > precOp(
//                    new AcaApproximateLuInverse<double>(discreteAcaLhs, delta));
//        // and wrap it in a PreconditionerBase object
//        // (the static cast is there because unspecifiedPrec() returns
//        // a ref-counted pointer to a subclass of PreconditionerBase
//        std::cout << "Created approximate inverse" << std::endl;
//        RCP<const Thyra::PreconditionerBase<double> > preconditioner =
//                Teuchos::rcp_static_cast<const Thyra::PreconditionerBase<double> >(
//                    Thyra::unspecifiedPrec(precOp));
//        // Now create a discrete linear operator with a solve() member function
//        invertibleDiscreteLhs = invertibleOpFactory.createOp();
//        Thyra::initializePreconditionedOp(
//                    invertibleOpFactory,
//                    trilinosDiscreteLhs,
//                    preconditioner,
//                    invertibleDiscreteLhs.ptr());
//        }
//        catch( const std::bad_cast &e ) {
//            std::cout << e.what() << std::endl;
//        }
