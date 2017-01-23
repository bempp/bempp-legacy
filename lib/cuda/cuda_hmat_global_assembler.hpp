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

#ifndef bempp_cuda_hmat_global_assembler_hpp
#define bempp_cuda_hmat_global_assembler_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/eigen_support.hpp"
#include "../common/types.hpp"

#include "../fiber/scalar_traits.hpp"
#include "../fiber/collection_of_kernels.hpp"

#include "../hmat/hmatrix.hpp"

#include <memory>
#include <vector>

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForIntegralOperators;
template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename GeometryFactor> class DefaultLocalAssemblerForIntegralOperatorsOnSurfaces;
template <typename BasisFunctionType> class Shapeset;
/** \endcond */

} // namespace Fiber

namespace Bempp {

/** \cond FORWARD_DECL */
template <typename ValueType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType> class Space;
template <typename BasisFunctionType, typename ResultType> class Context;
/** \endcond */

/** \ingroup weak_form_assembly_internal
 *  \brief CUDA HMAT-mode assembler.
 */
template <typename BasisFunctionType, typename ResultType>
class CudaHMatGlobalAssembler {
public:

  // TODO: How to determine the correct kernel type?
  typedef ResultType KernelType;

  typedef Fiber::Shapeset<BasisFunctionType> Shapeset;
  typedef DiscreteBoundaryOperator<ResultType> DiscreteBoundaryOp;
  typedef Fiber::LocalAssemblerForIntegralOperators<ResultType>
      LocalAssembler;

  static std::unique_ptr<DiscreteBoundaryOp> assembleDetachedWeakForm(
      const Space<BasisFunctionType> &testSpace,
      const Space<BasisFunctionType> &trialSpace,
      LocalAssembler &localAssembler,
      LocalAssembler &localAssemblerForAdmissibleBlocks,
      const Context<BasisFunctionType, ResultType> &context,
      int symmetry); // used to be "bool symmetric"; fortunately "true"
                     // is converted to 1 == SYMMETRIC

private:
  static void assembleDetachedWeakForm(
      const Space<BasisFunctionType> &testSpace,
      const Space<BasisFunctionType> &trialSpace,
      const shared_ptr<const Fiber::CollectionOfKernels<KernelType>> &kernel,
      const Shapeset &testShapeset, const Shapeset &trialShapeset,
      const std::vector<LocalAssembler*> &localAssemblers,
      const std::vector<LocalAssembler*> &localAssemblersForAdmissibleBlocks,
      const std::vector<const DiscreteBoundaryOp*> &sparseTermsToAdd,
      const std::vector<ResultType> &denseTermMultipliers,
      const std::vector<ResultType> &sparseTermMultipliers,
      const Context<BasisFunctionType, ResultType> &context, int symmetry,
      shared_ptr<hmat::DefaultHMatrixType<ResultType>> &hMatrix);
};

} // namespace Bempp

#endif
