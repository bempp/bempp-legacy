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
/** \endcond */
} // namespace Fiber

namespace Bempp {
/** \cond FORWARD_DECL */
template <typename ValueType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType> class Space;
template <typename BasisFunctionType, typename ResultType> class Context;

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>>
assembleDenseBlock(int rowStart, int rowEnd, int colStart, int colEnd,
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        Fiber::LocalAssemblerForIntegralOperators<ResultType>& assembler,
        const ParameterList& parameterList);




} // namespace Bempp

#endif
