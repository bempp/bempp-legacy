#ifndef bempp_context_hpp
#define bempp_context_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"
#include "../fiber/quadrature_strategy.hpp"
#include "assembly_options.hpp"
#include "discrete_boundary_operator_cache.hpp"

namespace Bempp
{

/** \cond FORWARD_DECL */
class GeometryFactory;
template <typename ValueType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class AbstractBoundaryOperator;
/** \endcond */

/** \ingroup weak_form_assembly
 *  \brief Assembly context.
 *
 *  This class manages the assembly of weak forms and evaluation of potentials.
 * 
 *  An assembly context consists of a quadrature strategy, which determines
 *  the way integrals are calculated, and assembly options, which control
 *  higher-level aspects of weak-form assembly, e.g. the use or not of
 *  acceleration algorithms such as ACA and the level of parallelism.
 *
 *  A call to the Context::getWeakForm() function returns the weak form of the
 *  abstract boundary operator passed in the argument to this function,
 *  assembled using the settings from the given Context. The Context may
 *  store the weak form in an internal cache; see the documentation of
 *  getWeakForm() for more details. */
template <typename BasisFunctionType, typename ResultType>
class Context
{
public:
    /** \brief Type of the appropriate instantiation of Fiber::QuadratureStrategy. */
    typedef Fiber::QuadratureStrategy<BasisFunctionType, ResultType, GeometryFactory>
    QuadratureStrategy;

    /** \brief Constructor.
     *
     *  \param[in] quadStrategy
     *    Quadrature strategy to be used for calculation of integrals occurring
     *    e.g. in the weak forms of boundary operators or in the definition of
     *    potential operators.
     *
     *  \param[in] assemblyOptions
     *    Further options influencing the weak-form assembly process. */
    Context(const shared_ptr<QuadratureStrategy>& quadStrategy,
            const AssemblyOptions& assemblyOptions);

    /** \brief Return the weak form of the specified abstract operator.
     *
     *  This function returns returns the weak form of the specified abstract
     *  boundary operator, calculated in accordance with the settings
     *  specified during the construction of the Context.
     *
     *  An important design principle in BEM++ is that abstract boundary
     *  operators are immutable after construction. The same is true for
     *  Context objects. Therefore, a Context stores newly calculated weak
     *  forms in an internal cache, unless a given abstract boundary operator
     *  does not provide a valid unique identifier (see
     *  AbstractBoundaryOperator::id()). As long as the weak form remains in
     *  cache, subsequent calls to Context::getWeakForm() with the same
     *  abstract boundary operator will not recalculate the weak form, but
     *  return the cached instance. It should, however, be noted that the cache
     *  does not maintain a persistent relationship with weak forms: it stores
     *  them as weak rather than shared pointers. Therefore, a weak form is
     *  deallocated as soon as the last reference to it *apart from that
     *  residing in the cache* goes out of scope. */
    shared_ptr<const DiscreteBoundaryOperator<ResultType> >
    getWeakForm(const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op) const {
        return m_cache.getWeakForm(*this, op);
    }

    /** \brief Return a reference to a copy of the AssemblyOptions object
     *  passed when constructing the Context. */
    const AssemblyOptions& assemblyOptions() const {
        return m_assemblyOptions;
    }

    /** \brief Return a reference to the QuadratureStrategy object
     *  passed when constructing the Context. */
    const QuadratureStrategy& quadStrategy() const {
        return *m_quadStrategy;
    }

private:
    shared_ptr<QuadratureStrategy> m_quadStrategy;
    AssemblyOptions m_assemblyOptions;

    mutable DiscreteBoundaryOperatorCache<BasisFunctionType, ResultType> m_cache;
};

} // namespace Bempp

#endif
