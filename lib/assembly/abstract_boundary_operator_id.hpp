#ifndef bempp_abstract_boundary_operator_id_hpp
#define bempp_abstract_boundary_operator_id_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"

#include <iostream>
#include <boost/functional/hash.hpp>

namespace Bempp
{

/** \ingroup abstract_boundary_operators
 *  \brief Base class for identifiers of an abstract boundary operator.
 *
 *  See AbstractBoundaryOperator::id() and Context::getWeakForm() for more
 *  information about the use of boundary-operator identifiers for caching
 *  weak forms.
 */
class AbstractBoundaryOperatorId
{
public:
    virtual ~AbstractBoundaryOperatorId() {}
    virtual size_t hash() const = 0;
    virtual void dump() const {}
    virtual bool isEqual(const AbstractBoundaryOperatorId& other) const = 0;

    bool operator==(const AbstractBoundaryOperatorId& other) const {
        return isEqual(other);
    }

    bool operator!=(const AbstractBoundaryOperatorId& other) const {
        return !operator==(other);
    }
};

} // namespace Bempp

namespace tbb
{

inline size_t tbb_hasher(
        const boost::shared_ptr<const Bempp::AbstractBoundaryOperatorId>& id)
{
    return id->hash();
}

inline size_t tbb_hasher(float id)
{
    return boost::hash<float>()(id);
}

inline size_t tbb_hasher(double id)
{
    return boost::hash<double>()(id);
}

} // namespace tbb

#include <tbb/concurrent_unordered_map.h>

namespace Bempp
{

// copied from Boost
template <typename T>
inline void tbb_hash_combine(size_t& seed, const T& v)
{
    seed ^= tbb::tbb_hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

}

#endif // bempp_abstract_boundary_operator_id_hpp
