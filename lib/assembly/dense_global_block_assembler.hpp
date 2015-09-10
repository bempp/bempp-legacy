#ifndef bempp_dense_global_block_assembler_hpp
#define bempp_dense_global_block_assembler_hpp

#include "../common/common.hpp"
#include "../common/eigen_support.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/scalar_traits.hpp"
#include "../common/types.hpp"

#include <memory>

namespace Fiber {
/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForIntegralOperators;
template <typename ResultType> class LocalAssemblerForPotentialOperators;
/** \endcond */
} // namespace Fiber

namespace Bempp {
/** \cond FORWARD_DECL */
class AssemblyOptions;
class EvaluationOptions;
template <typename ValueType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType> class Space;
template <typename BasisFunctionType, typename ResultType> class Context;
/** \endcond */
/** \ingroup weak_form_assembly_internal
 * * \brief Dense-mode assembler of matrix blocks.
 * */
template <typename BasisFunctionType, typename ResultType>
class DenseGlobalBlockAssembler {
public:
  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;
  typedef Fiber::LocalAssemblerForIntegralOperators<ResultType>
      LocalAssemblerForIntegralOperators;
  typedef Fiber::LocalAssemblerForPotentialOperators<ResultType>
      LocalAssemblerForPotentialOperators;
  static shared_ptr<DiscreteBoundaryOperator<ResultType>>
  assembleWeakForm(
      int colStart, int colEnd, int rowStart, int rowEnd,
      const Space<BasisFunctionType> &testSpace,
      const Space<BasisFunctionType> &trialSpace,
      LocalAssemblerForIntegralOperators &assembler,
      const ParameterList& parameterList);
};
} // namespace Bempp

#endif
