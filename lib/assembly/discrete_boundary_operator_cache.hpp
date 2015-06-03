#ifndef bempp_discrete_boundary_operator_cache_hpp
#define bempp_discrete_boundary_operator_cache_hpp

#include "../common/common.hpp"
#include "../common/deprecated.hpp"
#include "../common/shared_ptr.hpp"
#include <boost/scoped_ptr.hpp>
#include <vector>

namespace Bempp {

/** \cond FORWARD_DECL */
template <typename ResultType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType, typename ResultType>
class AbstractBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class Context;
/** \endcond */

/** \ingroup discrete_boundary_operators
 *  \brief Cache of discrete boundary operators.
 *
 *  \deprecated This class is no longer used and will be removed in a
 *  future version of BEM++.
 */
template <typename BasisFunctionType, typename ResultType>
class BEMPP_DEPRECATED DiscreteBoundaryOperatorCache {
public:
  /** \brief Constructor. */
  DiscreteBoundaryOperatorCache();

  /** \brief Destructor. */
  ~DiscreteBoundaryOperatorCache();

  /** \brief Return the weak form of the operator \p op.
   *
   *  This function first checks whether the weak form of \p op is already
   *  stored in the DiscreteBoundaryOperatorCache object. If it is, it is
   *  returned. Otherwise the weak form is assembled from scratch, possible
   *  stored in cache (if the operator \p op is cacheable) and returned to
   *  the caller. */
  BEMPP_DEPRECATED shared_ptr<const DiscreteBoundaryOperator<ResultType>>
  getWeakForm(
      const Context<BasisFunctionType, ResultType> &context,
      const AbstractBoundaryOperator<BasisFunctionType, ResultType> &op) const;

  /** \brief Return list of discrete operators currently stored in cache.
   *
   *  Can be used for debugging. */
  std::vector<shared_ptr<const DiscreteBoundaryOperator<ResultType>>>
  aliveOperators() const;

private:
  /** \cond PRIVATE */
  struct Impl;
  boost::scoped_ptr<Impl> m_impl;
  /** \endcond */
};

} // namespace Bempp

#endif
