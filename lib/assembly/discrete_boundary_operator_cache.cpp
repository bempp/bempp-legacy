#include "discrete_boundary_operator_cache.hpp"

#include "abstract_boundary_operator.hpp"
#include "abstract_boundary_operator_id.hpp"
#include "context.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include <boost/weak_ptr.hpp>
#include <tbb/concurrent_unordered_map.h>
#include <utility>

namespace Bempp
{

namespace
{

struct CompareSharedPtrsToConstAbstractBoundaryOperatorIds
{
    bool operator()(const shared_ptr<const Bempp::AbstractBoundaryOperatorId>& a,
                    const shared_ptr<const Bempp::AbstractBoundaryOperatorId>& b) const {
        return *a == *b;
    }
};

} // namespace

template <typename BasisFunctionType, typename ResultType>
struct DiscreteBoundaryOperatorCache<BasisFunctionType, ResultType>::Impl
{
    typedef tbb::concurrent_unordered_map<shared_ptr<const AbstractBoundaryOperatorId>,
    boost::weak_ptr<const DiscreteBoundaryOperator<ResultType> >,
    tbb::tbb_hash<shared_ptr<const AbstractBoundaryOperatorId> >,
    CompareSharedPtrsToConstAbstractBoundaryOperatorIds >
    DiscreteBoundaryOperatorMap;

    mutable DiscreteBoundaryOperatorMap discreteOps;
};

template <typename BasisFunctionType, typename ResultType>
DiscreteBoundaryOperatorCache<BasisFunctionType, ResultType>::
DiscreteBoundaryOperatorCache() :
  m_impl(new Impl)
{
}

template <typename BasisFunctionType, typename ResultType>
DiscreteBoundaryOperatorCache<BasisFunctionType, ResultType>::
~DiscreteBoundaryOperatorCache()
{
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType> >
DiscreteBoundaryOperatorCache<BasisFunctionType, ResultType>::getWeakForm(
        const Context<BasisFunctionType, ResultType>& context,
        const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op) const
{
    // std::cout << "Weak form of operator of type " << typeid(op).name()
    //           << " requested from cache" << std::endl;
    typedef DiscreteBoundaryOperator<ResultType> DiscreteOp;

    shared_ptr<const AbstractBoundaryOperatorId> id = op.id();

    if (!id) { // operator is not cacheable
        // std::cout << "Noncacheable discrete operator requested" << std::endl;
        return op.assembleWeakForm(context);
    }

    // Check if a weak pointer to discrete operator exists in cache
    typename Impl::DiscreteBoundaryOperatorMap::const_iterator it =
            m_impl->discreteOps.find(id);
    if (it != m_impl->discreteOps.end()) {
        // std::cout << "Weak pointer to discrete operator found in cache" << std::endl;
        // Check if the weak pointer still refers to an existing object
        if (shared_ptr<const DiscreteOp> cachedDiscreteOp = it->second.lock()) {
            // std::cout << "Weak pointer is valid" << std::endl;
            return cachedDiscreteOp;
        } else {
            // std::cout << "Weak pointer is invalid, reconstructing discrete operator" << std::endl;
            shared_ptr<const DiscreteOp> discreteOp = op.assembleWeakForm(context);
            // TODO: THIS IS NOT NOT THREAD-SAFE!!!
            it->second = discreteOp;
            return it->second.lock();
        }
    }

    // std::cout << "Discrete operator not found in cache -- constructing from scratch"
    //           << std::endl;
    // Discrete operator not found in cache, construct it...
    shared_ptr<const DiscreteOp> discreteOp = op.assembleWeakForm(context);
    // and attempt to insert it into the map
    std::pair<typename Impl::DiscreteBoundaryOperatorMap::iterator, bool>
            result = m_impl->discreteOps.insert(std::make_pair(id, discreteOp));
    // Return pointer to the discrete operator that ended up in the map.
    return result.first->second.lock();

    // (Note: if another thread was faster in trying to insert a discrete
    // operator into the map, the operator we constructed (pointed by
    // discreteOp) will be destructed now, as discreteOp goes out of scope.)
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(DiscreteBoundaryOperatorCache);

} // namespace Bempp
