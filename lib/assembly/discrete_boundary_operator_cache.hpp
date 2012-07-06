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
    bool operator()(const boost::shared_ptr<const Bempp::AbstractBoundaryOperatorId>& a,
                    const boost::shared_ptr<const Bempp::AbstractBoundaryOperatorId>& b) const {
        return *a == *b;
    }
};

template <typename BasisFunctionType, typename ResultType>
class DiscreteBoundaryOperatorCache
{
public:
    DiscreteBoundaryOperatorCache() {}

    shared_ptr<const DiscreteBoundaryOperator<ResultType> >
    getWeakForm(const Context<BasisFunctionType, ResultType>& context,
                const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op) const;

private:
    typedef tbb::concurrent_unordered_map<shared_ptr<const AbstractBoundaryOperatorId>,
    boost::weak_ptr<const DiscreteBoundaryOperator<ResultType> >,
    tbb::tbb_hash<shared_ptr<const AbstractBoundaryOperatorId> >,
    CompareSharedPtrsToConstAbstractBoundaryOperatorIds >
    DiscreteBoundaryOperatorMap;

    mutable DiscreteBoundaryOperatorMap m_discreteOps;
};

} // namespace Bempp

#endif
