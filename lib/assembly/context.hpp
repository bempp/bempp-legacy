#ifndef bempp_context_hpp
#define bempp_context_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"
#include "../fiber/quadrature_strategy.hpp"
#include "../common/global_parameters.hpp"
#include "../common/types.hpp"
#include "assembly_options.hpp"
#include "../cuda/cuda_options.hpp"
#include "discrete_boundary_operator_cache.hpp"

namespace Bempp {

/** \cond FORWARD_DECL */
class GeometryFactory;
template <typename ValueType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType, typename ResultType>
class AbstractBoundaryOperator;
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
 */
template <typename BasisFunctionType, typename ResultType> class Context {
public:
  /** \brief Type of the appropriate instantiation of Fiber::QuadratureStrategy.
   */
  typedef Fiber::QuadratureStrategy<BasisFunctionType, ResultType,
                                    GeometryFactory> QuadratureStrategy;

  /** \brief Constructor.
   *
   *  \param[in] quadStrategy
   *    Quadrature strategy to be used for calculation of integrals occurring
   *    e.g. in the weak forms of boundary operators or in the definition of
   *    potential operators.
   *
   *  \param[in] assemblyOptions
   *    Further options influencing the weak-form assembly process.
   *
   *  \param[in] globalParameterList
   *    New parameter structure that will supersede the other option objects. */
  Context(const shared_ptr<const QuadratureStrategy> &quadStrategy,
          const AssemblyOptions &assemblyOptions,
          const ParameterList &globalParameterList =
              GlobalParameters::parameterList());

  /** \brief Constructor */
  Context();

  /** \brief Constructor.
   *
   *  \param[in] globalParameterList
   *     Parameter list that contains the BEM++ options. */
  explicit Context(const ParameterList &globalParameterList);

  /** \brief Return the discrete weak form of the specified abstract operator.
   *
   *  This function returns returns the weak form of the specified abstract
   *  boundary operator, calculated in accordance with the settings
   *  specified during the construction of the Context.
   *
   *  \deprecated This function is deprecated and will be removed in a future
   *  version of BEM++. To get the discrete weak form of an
   *  AbstractBoundaryOperator, use its <tt>assembleWeakForm()</tt> method.
   *  Alternatively, call the <tt>weakForm()</tt> method of a
   *  BoundaryOperator object wrapping the AbstractBoundaryOperator in
   *  question. */
  BEMPP_DEPRECATED
  shared_ptr<const DiscreteBoundaryOperator<ResultType>> getWeakForm(
      const AbstractBoundaryOperator<BasisFunctionType, ResultType> &op) const;

  /** \brief Return a reference to a copy of the AssemblyOptions object
   *  passed when constructing the Context. */
  const AssemblyOptions &assemblyOptions() const { return m_assemblyOptions; }

  /** \brief Return a reference to the QuadratureStrategy object
   *  passed when constructing the Context. */
  shared_ptr<const QuadratureStrategy> quadStrategy() const {
    return m_quadStrategy;
  }

  /** \brief Return a reference to a copy of the CudaOptions object */
  const CudaOptions &cudaOptions() const { return m_cudaOptions; }

  /** \brief Const version of \p globalParameterList. */

  const ParameterList &globalParameterList() const {
    return m_globalParameterList;
  }

private:
  shared_ptr<const QuadratureStrategy> m_quadStrategy;
  AssemblyOptions m_assemblyOptions;
  CudaOptions m_cudaOptions;
  ParameterList m_globalParameterList;
};

} // namespace Bempp

#endif
