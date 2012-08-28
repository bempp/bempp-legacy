#ifndef bempp_discrete_boundary_operator_cache_hpp
#define bempp_discrete_boundary_operator_cache_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"

#include "abstract_boundary_operator_id.hpp"

#include <boost/weak_ptr.hpp>
#include <tbb/concurrent_unordered_map.h>

namespace Bempp
{

class AbstractBoundaryOperatorId;
template <typename ResultType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class AbstractBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class Context;

struct CompareSharedPtrsToConstAbstractBoundaryOperatorIds
{
    bool operator()(const shared_ptr<const Bempp::AbstractBoundaryOperatorId>& a,
                    const shared_ptr<const Bempp::AbstractBoundaryOperatorId>& b) const {
        return *a == *b;
    }
};

/** \ingroup discrete_boundary_operators
 *  \brief Cache of discrete boundary operators.
 */
template <typename BasisFunctionType, typename ResultType>
class DiscreteBoundaryOperatorCache
{
public:
    /** \brief Constructor. */
    DiscreteBoundaryOperatorCache() {}

    /** \brief Return the weak form of the operator \p op.
     *
     *  This function first checks whether the weak form of \p op is already
     *  stored in the DiscreteBoundaryOperatorCache object. If it is, it is
     *  returned. Otherwise the weak form is assembled from scratch, possible
     *  stored in cache (if the operator \p op is cacheable) and returned to
     *  the caller. */
    shared_ptr<const DiscreteBoundaryOperator<ResultType> >
    getWeakForm(const Context<BasisFunctionType, ResultType>& context,
                const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op) const;

private:
    /** \cond PRIVATE */
    typedef tbb::concurrent_unordered_map<shared_ptr<const AbstractBoundaryOperatorId>,
    boost::weak_ptr<const DiscreteBoundaryOperator<ResultType> >,
    tbb::tbb_hash<shared_ptr<const AbstractBoundaryOperatorId> >,
    CompareSharedPtrsToConstAbstractBoundaryOperatorIds >
    DiscreteBoundaryOperatorMap;

    mutable DiscreteBoundaryOperatorMap m_discreteOps;
    /** \endcond */
};

} // namespace Bempp

#endif
