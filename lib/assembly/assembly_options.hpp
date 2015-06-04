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

#ifndef bempp_assembly_options_hpp
#define bempp_assembly_options_hpp

#include "../common/common.hpp"

#include "aca_options.hpp"

#include "../common/deprecated.hpp"
#include "../fiber/opencl_options.hpp"
#include "../fiber/parallelization_options.hpp"
#include "../fiber/verbosity_level.hpp"

namespace Bempp {

using Fiber::OpenClOptions;
using Fiber::ParallelizationOptions;
using Fiber::VerbosityLevel;

/** \ingroup weak_form_assembly
 *  \brief Options determining how weak-form assembly is done.
 */
class AssemblyOptions {
public:
  enum Value { AUTO = -1, NO = 0, YES = 1 };

  AssemblyOptions();

  /** @name Assembly mode
    @{ */

  /** \brief Possible assembly modes for weak forms of boundary integral
   * operators. */
  enum Mode {
    /** \brief Assemble dense matrices. */
    DENSE,
    /** \brief Assemble hierarchical matrices using adaptive cross approximation
       (ACA). */
    ACA,
    /** \brief Assemble hierarchical matrices using the HMat library. */
    HMAT
  };

  /** \brief Use dense-matrix representations of weak forms of boundary integral
   *operators.
   *
   *  This is the default assembly mode. */
  void switchToDenseMode();

  /** \brief Use adaptive cross approximation (ACA) to obtain
   *hierarchical-matrix
   *  representations of weak forms of boundary integral operators.
   *
   *  \param[in] acaOptions Parameters influencing the ACA algorithm. */
  void switchToAcaMode(const AcaOptions &acaOptions);

  /** \brief Assemble using the HMat hierarchical matrix library. */
  void switchToHMatMode();

  /** \brief Use dense-matrix representations of weak forms of boundary integral
   *operators.
   *
   *  \deprecated Use switchToDenseMode() instead. */
  BEMPP_DEPRECATED void switchToDense();

  /** \brief Use adaptive cross approximation (ACA) to obtain
   *hierarchical-matrix
   *  representations of weak forms of boundary integral operators.
   *
   *  \deprecated Use switchToAcaMode() instead. */
  BEMPP_DEPRECATED void switchToAca(const AcaOptions &acaOptions);

  /** \brief Current assembly mode.
   *
   *  The assembly mode can be changed by calling switchToDenseMode() or
   *  switchToAcaMode(). */
  Mode assemblyMode() const;

  /** \brief Return the current adaptive cross approximation (ACA) settings.
   *
   *  \note These settings are only used in the ACA assembly mode, i.e. when
   *  assemblyMode() returns ACA. */
  const AcaOptions &acaOptions() const;

  /** @}
    @name Parallelization
    @{ */

  // Temporarily removed (OpenCl support is broken).
  // void enableOpenCl(const OpenClOptions& openClOptions);
  // void disableOpenCl();

  /** \brief Set the maximum number of threads used during the assembly.
   *
   *  \p maxThreadCount must be a positive number or \p AUTO. In the latter
   *  case the number of threads is determined automatically. */
  void setMaxThreadCount(int maxThreadCount);

  /** \brief Set the maximum number of threads used during the assembly.
   *
   *  \deprecated Use setMaxThreadCount() instead. */
  BEMPP_DEPRECATED void switchToTbb(int maxThreadCount = AUTO);

  /** \brief Return current parallelization options. */
  const ParallelizationOptions &parallelizationOptions() const;

  /** @}
    @name Verbosity
    */

  /** \brief Set the verbosity level.
   *
   *  This setting determines the amount of information printed out by
   *  functions from BEM++. */
  void setVerbosityLevel(VerbosityLevel::Level level);

  /** \brief Return the verbosity level. */
  VerbosityLevel::Level verbosityLevel() const;

  /** @}
    @name Miscellaneous
    @{ */

  /** \brief Specify whether singular integrals are cached during weak-form
   *assembly.
   *
   *  If <tt>value == true</tt>, singular integrals are precalculated
   *  and stored in a cache before filling the matrix of the discretized weak
   *  form of a singular boundary operator. Otherwise these integrals
   *  are evaluated as needed during the assembly of the weak form.
   *
   *  By default, singular integral caching is enabled. */
  void enableSingularIntegralCaching(bool value = true);

  /** \brief Return whether singular integrals should be cached during weak-form
   *assembly.
   *
   *  See enableSingularIntegralCaching() for more information. */
  bool isSingularIntegralCachingEnabled() const;

  /** \brief Specify whether discrete weak forms of local operators should be
   *  stored in sparse format.
   *
   *  If <tt>value == true</tt> (default), discretized local operators are
   *  stored as sparse matrices, otherwise as dense matrices. */
  void enableSparseStorageOfLocalOperators(bool value = true);

  /** \brief Return whether discrete weak forms of local operators should be
   *  stored in sparse format.
   *
   *  See enableSparseStorageOfLocalOperators() for more information. */
  bool isSparseStorageOfLocalOperatorsEnabled() const;

  /** \brief Specify whether discrete weak forms of local operators should be
   *  stored in sparse format.
   *
   *  \deprecated This function is deprecated and will be removed in a future
   *  version of BEM++. Use enableSparseStorageOfLocalOperators() instead.
   */
  void BEMPP_DEPRECATED enableSparseStorageOfMassMatrices(bool value = true);

  /** \brief Return whether mass matrices should be stored in sparse format.
   *
   *  See enableSparseStorageOfMassMatrices() for more information. */
  bool BEMPP_DEPRECATED isSparseStorageOfMassMatricesEnabled() const;

  /** \brief Enable or disable joint assembly of integral-operator
   *superpositions.
   *
   *  If <tt>value == true</tt>, the discrete weak forms of superpositions of
   *  elementary integral operators, such as 5. * SLP + 2. * DLP, will
   *  assembled together (for example, in a single run of ACA) and stored as a
   *  single H-matrix or dense matrix. By default joint assembly is disabled,
   *  which means that the discrete weak form of each elementary operator is
   *  stored separately. */
  void enableJointAssembly(bool value = true);

  /** \brief Return whether joint assembly of integral-operator superpositions
   *  is enabled.
   *
   * See enableJointAssembly() for more information. */
  bool isJointAssemblyEnabled() const;

  /** \brief Specify whether BLAS matrix multiplication routines should be
   *  used during evaluation of elementary integrals.
   *
   *  If this option is set to AUTO (default), BLAS is used in the evaluation
   *  of integrals occurring in the weak forms of operators whose test or
   *  trial space is composed of quadratic or higher-order elements; for the
   *  remaining operators, BLAS routines are not used. You can force
   *  BLAS-based integration routines to be used always (or never) by setting
   *  this option to \c YES (or \c NO).
   */
  void enableBlasInQuadrature(Value value = AUTO);

  /** \brief Indicate whether BLAS matrix multiplication routines are used
   *  during evaluation of elementary integrals.
   *
   *  See enableBlasInQuadrature() for more information. */
  Value isBlasEnabledInQuadrature() const;

  /** \brief Instruct the ACA assembler to use the same quadrature order for
   *  the evaluation of all regular integrals over pairs of test and trial
   *  elements belonging to a single block cluster.
   *
   *  Use of different quadrature rules to evaluate regular integrals making
   *  up the entries of a single block of an H-matrix may lead to an increase
   *  in the block rank, and hence memory consumption, if these quadrature
   *  rules are not accurate enough (if their error is greater than the ACA
   *  tolerance \f$\epsilon\f$). However, evaluation of the entries of large
   *  blocks using a uniform quadrature rule may be unnecessarily
   *  time-consuming.
   *
   *  If this option is set to \c false and one uses a quadrature strategy
   *  that takes the interlement distance into account when selecting
   *  quadrature order, the true interelement distance is calculated and used
   *  for each pair of elements. Otherwise, instead of the distance between
   *  elements the quadrature-rule selector is passed the distance between
   *  the clusters to which these elements belong.
   *
   *  By default, this option is set to \c true.
   *  \note In HMat always true. */
  void makeQuadratureOrderUniformInEachCluster(bool value = true);

  /** \brief Return whether the same quadrature order is used for all regular
   *  integrals over pairs of test and trial elements belonging to a single
   *  block cluster.
   *
   *  See makeQuadratureOrderUniformInEachCluster() for more information. */
  bool isQuadratureOrderUniformInEachCluster() const;

  /** @} */

private:
  /** \cond */
  Mode m_assemblyMode;
  AcaOptions m_acaOptions;
  ParallelizationOptions m_parallelizationOptions;
  VerbosityLevel::Level m_verbosityLevel;
  bool m_singularIntegralCaching;
  bool m_sparseStorageOfLocalOperators;
  bool m_jointAssembly;
  bool m_uniformQuadrature;
  Value m_blasInQuadrature;
  /** \endcond */
};

} // namespace Bempp

#endif
