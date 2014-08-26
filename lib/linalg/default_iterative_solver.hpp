// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef bempp_default_iterative_solver_hpp
#define bempp_default_iterative_solver_hpp

#include "../common/common.hpp"
#include "bempp/common/config_trilinos.hpp"

#ifdef WITH_TRILINOS

#include "solver.hpp"

#include "belos_solver_wrapper_fwd.hpp" // for default parameter lists
#include "preconditioner.hpp"

#include "../common/deprecated.hpp"

#include <boost/scoped_ptr.hpp>

namespace Thyra {
/** \cond FORWARD_DECL */
template <typename ValueType> class PreconditionerBase;
/** \endcond */
} // namespace Thyra

namespace Bempp {

/** \ingroup linalg
  * \brief Default iterative solver for boundary integral equations.
  *
  * This class provides an interface to various iterative solvers available via
  * the Stratimikos interface to the Belos solver family from Trilinos (see <a
  * href="http://trilinos.sandia.gov/packages/docs/r10.10/packages/stratimikos/doc/html/index.html">Stratimikos
  * documentation</a>). Convergence can be tested either in range space or in
  * the dual space to the range space. A standard Galerkin discretisation of
  * the form \f$Ax=b\f$, maps into the dual space of the range of the operator.
  * By choosing to test in the range space the equation \f$M^\dagger
  * Ax=M^\dagger b\f$ is solved, where \f$M\f$ is the mass matrix, mapping from
  * the range space into its dual and \f$M^\dagger\f$ is its pseudoinverse.
  *
  */
template <typename BasisFunctionType, typename ResultType>
class DefaultIterativeSolver : public Solver<BasisFunctionType, ResultType> {
public:
  typedef Solver<BasisFunctionType, ResultType> Base;

  /** \brief Constructor of the <tt>DefaultIterativeSolver</tt> class.
    *
    * \param[in] boundaryOp
    *   Non-blocked boundary operator.
    * \param[in] mode
    *   Convergence test mode. Default:
    *<tt>TEST_CONVERGENCE_IN_DUAL_TO_RANGE</tt>.
    *
    */
  DefaultIterativeSolver(
      const BoundaryOperator<BasisFunctionType, ResultType> &boundaryOp,
      ConvergenceTestMode::Mode mode =
          ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);

  /** \brief Constructor of the <tt>DefaultIterativeSolver</tt> class.
    *
    * \param[in] boundaryOp
    *   Blocked boundary operator.
    * \param[in] mode
    *   Convergence test mode. Default:
    *<tt>TEST_CONVERGENCE_IN_DUAL_TO_RANGE</tt>.
    *
    */
  DefaultIterativeSolver(
      const BlockedBoundaryOperator<BasisFunctionType, ResultType> &boundaryOp,
      ConvergenceTestMode::Mode mode =
          ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);

  virtual ~DefaultIterativeSolver();

  /** \brief Define a preconditioner.
    *
    * The preconditioner is passed on to the Belos solver.
    *
    * \param[in] preconditioner
    *
    * \note This function has no effect if it is called after
    *initializeSolver().
    *
    * \deprecated Do not use this function in new code. Instead, pass the
    *   preconditioner to initializeSolver().
    */
  BEMPP_DEPRECATED void
  setPreconditioner(const Preconditioner<ResultType> &preconditioner);

  /** \brief Initialize a Belos iterative solver.
    *
    * \param[in] paramList
    *   Parameter lists can be read in from XML files or defined in code. Use
    *   defaultGmresParameterList() and defaultCgParameterList() to construct
    *   default parameter lists for the GMRES and CG solvers.
    */
  void initializeSolver(const Teuchos::RCP<Teuchos::ParameterList> &paramList);

  /** \brief Initialize a preconditioned Belos iterative solver.
    *
    * \param[in] paramList
    *   Parameter lists can be read in from XML files or defined in code. Use
    *   defaultGmresParameterList() and defaultCgParameterList() to construct
    *   default parameter lists for the GMRES and CG solvers.
    * \param[in] preconditioner
    *   Preconditioner to be used by the solver.
    */
  void initializeSolver(const Teuchos::RCP<Teuchos::ParameterList> &paramList,
                        const Preconditioner<ResultType> &preconditioner);

private:
  virtual Solution<BasisFunctionType, ResultType> solveImplNonblocked(
      const GridFunction<BasisFunctionType, ResultType> &rhs) const;
  virtual BlockedSolution<BasisFunctionType, ResultType> solveImplBlocked(
      const std::vector<GridFunction<BasisFunctionType, ResultType>> &rhs)
      const;

private:
  struct Impl;
  boost::scoped_ptr<Impl> m_impl;
};

} // namespace Bempp

#endif // WITH_TRILINOS

#endif
