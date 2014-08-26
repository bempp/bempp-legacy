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

#ifndef bempp_solution_base_hpp
#define bempp_solution_base_hpp

#include "bempp/common/config_trilinos.hpp"
#include "../common/common.hpp"

#include "../common/scalar_traits.hpp"

// TODO: support lack of trilinos
#ifdef WITH_TRILINOS
#include <Thyra_SolveSupportTypes.hpp>
#endif // WITH_TRILINOS

#include <string>

namespace Bempp {

struct SolutionStatus {
  enum Status {
    CONVERGED,
    UNCONVERGED,
    UNKNOWN
  };
};

/** \ingroup linalg
  * \brief The base class for the Solution and BlockedSolution container
  *classes.
  *
  * This class holds various information about the solution of a BEM
  *computation. It derives from
  * <a
  *href="http://trilinos.sandia.gov/packages/docs/r10.10/packages/thyra/doc/html/structThyra_1_1SolveStatus.html">Thyra::SolveStatus</a>
  * if Trilinos is enabled.
  *
  * <tt>MagnitudeType</tt> is identical to <tt>float</tt> if
  * <tt>ResultType</tt> is <tt>float</tt> or <tt>complex<float></tt>. It is
  * double if <tt>ResultType</tt> is <tt>double</tt> or
  * <tt>complex<double></tt>
  */
template <typename BasisFunctionType, typename ResultType> class SolutionBase {
public:
  typedef typename ScalarTraits<ResultType>::RealType MagnitudeType;

#ifdef WITH_TRILINOS
  /** \brief Constructor */
  explicit SolutionBase(const Thyra::SolveStatus<MagnitudeType> status);
#endif // WITH_TRILINOS
  /** \brief Constructor */
  explicit SolutionBase(SolutionStatus::Status status,
                        MagnitudeType achievedTolerance = unknownTolerance(),
                        std::string message = "");

  static MagnitudeType unknownTolerance() { return MagnitudeType(-1.); }

  /** \brief Return status of the linear solve. */
  SolutionStatus::Status status() const;

  /** \brief Maximum final tolerance achieved by the linear solve.
   *
   *  A value of unknownTolerance() means that even an estimate of
   *  the the final value of the tolerance is unknown.  This is the
   *  typical value returned by direct solvers. */
  MagnitudeType achievedTolerance() const;

  /** \brief Iteration count if iterative solver is used.
   *  If a direct solver is used, this defaults to -1.
   */
  int iterationCount() const;

  /** \brief Message returned by the solver. */
  std::string solverMessage() const;

#ifdef WITH_TRILINOS
  /** \brief Extra status parameter returned by the solver.
   *
   *  \note Contents of this list are solver-dependent. */
  Teuchos::RCP<Teuchos::ParameterList> extraParameters() const;
#endif // WITH_TRILINOS

private:
  SolutionStatus::Status m_status;
  MagnitudeType m_achievedTolerance;
  std::string m_message;
  int m_iterationCount;
#ifdef WITH_TRILINOS
  Teuchos::RCP<Teuchos::ParameterList> m_extraParameters;
#endif // WITH_TRILINOS
};

} // namespace Bempp

#endif
