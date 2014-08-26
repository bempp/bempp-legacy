#ifndef bempp_dense_global_assembler_hpp
#define bempp_dense_global_assembler_hpp

#include "../common/common.hpp"

#include <memory>

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForIntegralOperators;
/** \endcond */

} // namespace Fiber

namespace Bempp {

/** \cond FORWARD_DECL */
class AssemblyOptions;
template <typename ValueType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType> class Space;
template <typename BasisFunctionType, typename ResultType> class Context;
/** \endcond */

/** \ingroup weak_form_assembly_internal
 *  \brief Dense-mode assembler.
 */
template <typename BasisFunctionType, typename ResultType>
class DenseGlobalAssembler {
public:
  typedef Fiber::LocalAssemblerForIntegralOperators<ResultType>
  LocalAssemblerForIntegralOperators;

  static std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
  assembleDetachedWeakForm(
      const Space<BasisFunctionType> &testSpace,
      const Space<BasisFunctionType> &trialSpace,
      LocalAssemblerForIntegralOperators &assembler,
      const Context<BasisFunctionType, ResultType> &context);
};

} // namespace Bempp

#endif
