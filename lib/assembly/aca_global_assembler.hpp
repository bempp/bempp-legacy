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


#ifndef bempp_aca_global_assembler_hpp
#define bempp_aca_global_assembler_hpp

#include "../common/common.hpp"

#include "../common/armadillo_fwd.hpp"
#include "../common/shared_ptr.hpp"
#include "../fiber/scalar_traits.hpp"

#include <memory>
#include <vector>

class Epetra_CrsMatrix;

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForOperators;
template <typename ResultType> class LocalAssemblerForPotentialOperators;
/** \endcond */

} // namespace Fiber

namespace Bempp
{

/** \cond FORWARD_DECL */
class AssemblyOptions;
class EvaluationOptions;
template <typename ValueType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType> class Space;
template <typename BasisFunctionType, typename ResultType> class Context;
/** \endcond */

/** \ingroup weak_form_assembly_internal
 *  \brief ACA-mode assembler.
 */
template <typename BasisFunctionType, typename ResultType>
class AcaGlobalAssembler
{
    typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;

public:
    typedef DiscreteBoundaryOperator<ResultType> DiscreteBndOp;
    typedef Fiber::LocalAssemblerForOperators<ResultType>
    LocalAssemblerForBoundaryOperators;
    typedef LocalAssemblerForBoundaryOperators LocalAssembler; // deprecated
    typedef Fiber::LocalAssemblerForPotentialOperators<ResultType>
    LocalAssemblerForPotentialOperators;

    static std::auto_ptr<DiscreteBndOp> assembleDetachedWeakForm(
            const Space<BasisFunctionType>& testSpace,
            const Space<BasisFunctionType>& trialSpace,
            const std::vector<LocalAssemblerForBoundaryOperators*>& localAssemblers,
            const std::vector<const DiscreteBndOp*>& sparseTermsToAdd,
            const std::vector<ResultType>& denseTermMultipliers,
            const std::vector<ResultType>& sparseTermMultipliers,
            const Context<BasisFunctionType, ResultType>& context,
            int symmetry);

    static std::auto_ptr<DiscreteBndOp> assembleDetachedWeakForm(
            const Space<BasisFunctionType>& testSpace,
            const Space<BasisFunctionType>& trialSpace,
            LocalAssemblerForBoundaryOperators& localAssembler,
            const Context<BasisFunctionType, ResultType>& context,
            int symmetry); // used to be "bool symmetric"; fortunately "true"
                           // is converted to 1 == SYMMETRIC

    static std::auto_ptr<DiscreteBndOp> assemblePotentialOperator(
            const arma::Mat<CoordinateType>& points,
            const Space<BasisFunctionType>& trialSpace,
            const std::vector<LocalAssemblerForPotentialOperators*>& localAssemblers,
            const std::vector<ResultType>& termMultipliers,
            const EvaluationOptions& options);

    static std::auto_ptr<DiscreteBndOp> assemblePotentialOperator(
            const arma::Mat<CoordinateType>& points,
            const Space<BasisFunctionType>& trialSpace,
            LocalAssemblerForPotentialOperators& localAssembler,
            const EvaluationOptions& options);
};

} // namespace Bempp

#endif
