#ifndef bempp_dense_global_assembler_hpp
#define bempp_dense_global_assembler_hpp

#include "../common/common.hpp"
#include "../common/armadillo_fwd.hpp"
#include "../common/scalar_traits.hpp"

#include <memory>

namespace Fiber
{
    /** \cond FORWARD_DECL */
    template <typename ResultType> class LocalAssemblerForIntegralOperators;
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
     * * \brief Dense-mode assembler.
     * */
    template <typename BasisFunctionType, typename ResultType>
        class DenseGlobalAssembler
        {
            public:
                typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;
                typedef Fiber::LocalAssemblerForIntegralOperators<ResultType>
                    LocalAssemblerForIntegralOperators;
                typedef Fiber::LocalAssemblerForPotentialOperators<ResultType>
                    LocalAssemblerForPotentialOperators;
                static std::unique_ptr<DiscreteBoundaryOperator<ResultType> >
                    assembleDetachedWeakForm(
                            const Space<BasisFunctionType>& testSpace,
                            const Space<BasisFunctionType>& trialSpace,
                            LocalAssemblerForIntegralOperators& assembler,
                            const Context<BasisFunctionType, ResultType>& context);
                static std::unique_ptr<DiscreteBoundaryOperator<ResultType> >
                    assemblePotentialOperator(
                            const arma::Mat<CoordinateType>& points,
                            const Space<BasisFunctionType>& trialSpace,
                            LocalAssemblerForPotentialOperators& assembler,
                            const EvaluationOptions& options);
        };
} // namespace Bempp

#endif
