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

#ifndef bempp_aca_options_hpp
#define bempp_aca_options_hpp

#include "../common/common.hpp"

#include <string>

namespace Bempp {

/** \ingroup weak_form_assembly
 *  \brief Adaptive cross approximation (ACA) parameters.
 */
class AcaOptions {
public:
  /** \brief Initialize ACA parameters to default values. */
  AcaOptions();

  /** \brief Estimate of the desired approximation accuracy.
   *
   *  Default value: 1e-4. */
  double eps;

  /** \brief Cluster-pair admissibility parameter.
   *
   *  Default value: 1.2. */
  double eta;

  /** \brief Minimum size of blocks approximated with ACA.
   *
   *  Matrix blocks whose smaller side is less than this value will be stored
   *  in the dense format (ACA will not be attempted on these blocks).
   *
   *  Default value: 16. */
  unsigned int minimumBlockSize;

  /** \brief Maximum allowed block size.
   *
   *  Matrix blocks with side larger than this value will be split. This can
   *  be used to limit the time devoted to (serial) processing of an
   *  individual block, especially when a large number of threads is used.
   *
   *  Default value: UINT_MAX. */
  unsigned int maximumBlockSize;

  /** \brief Maximum rank of blocks stored in the low-rank format.
   *
   *  Blocks judged to have higher rank will be stored in the dense format.
   *
   *  Default value: UINT_MAX. */
  unsigned int maximumRank;

  /** \brief ACA assembly mode. See documentation of the member \p mode for
   *  more information. */
  enum AcaAssemblyMode {
    MIN_ASSEMBLY_MODE,
    GLOBAL_ASSEMBLY = MIN_ASSEMBLY_MODE,
    LOCAL_ASSEMBLY,
    HYBRID_ASSEMBLY,
    MAX_ASSEMBLY_MODE = HYBRID_ASSEMBLY
  };

  /** \brief ACA assembly mode.
   *
   *  This parameter lets you choose between three ways of discretizing
   *  operators using ACA. Suppose you have an integral operator \f$\mathcal
   *  A\f$ that you want to discretize using a test space \f$U\f$ and and
   *  trial space \f$V\f$.
   *
   *  1. If \p mode is set to \p GLOBAL_ASSEMBLY (default), ACA is used in
   *     the most straightforward way possible, i.e. to approximate blocks of
   *     the operator \f$\mathcal A\f$ discretized with the basis functions of
   *     \f$U\f$ and \f$V\f$.
   *
   *  2. If \p mode is set to \p LOCAL_ASSEMBLY, the discretization \f$A\f$ of
   *     \f$\mathcal A\f$ is constructed as a sum
   *
   *     \f[
   *         A = \sum_{i=1}^n \alpha_i P_i A_i Q_i,
   *     \f]
   *
   *     where \f$n\f$ is a small integer, \f$\alpha_i\f$ scalar
   *     coefficients, \f$P_i\f$ and \f$Q_i\f$ sparse matrices and \f$A_i\f$
   *     H-matrices representing integral operators discretized with basis
   *     functions whose support extends over single elements only (for
   *     example functions linear on a single element and zero everywhere
   *     else). In most practical cases all matrices \f$A_i\f$ are equal (and
   *     obviously only a single copy is stored in memory). This assembly
   *     mode takes advantage of the fact that ACA becomes more efficient as
   *     the support of test and trial functions is reduced. Thus, assembly
   *     in LOCAL_ASSEMBLY mode is often much faster than in GLOBAL_ASSEMBLY
   *     mode. However, this comes at the price of increased memory use,
   *     since the matrices \f$A_i\f$ are typically much bigger than \f$A\f$.
   *     In addition, the
   *     DiscreteBoundaryOperator::asDiscreteAcaBoundaryOperator() function
   *     is not currently implemented for products of operators, so discrete
   *     operators created in the \p LOCAL_ASSEMBLY modes cannot be converted
   *     into single H-matrices. Thus, they cannot be used for H-LU
   *     preconditioning.
   *
   *     More information about the matrices \f$P_i\f$ and \f$Q_i\f$ can be
   *     found in the documentation of the SyntheticIntegralOperator class.
   *     See also the documentation of non-member constructors of particular
   *     integral operators, such as laplace3dSingleLayerBoundaryOperator()
   *     and helmholtz3dHypersingularBoundaryOperator(), for explicit
   *     expressions used to represent these operators in the \p
   *     LOCAL_ASSEMBLY mode.
   *
   *  3. The \p HYBRID_ASSEMBLY mode is available for operators whose
   *     weak form can be written as
   *
   *     \f[
   *         \int_\Gamma \int_\Sigma
   *         f(x) \, K(x, y) \, g(y) \,
   *         \mathrm{d}\Gamma(x) \,\mathrm{d}\Sigma(y),
   *     \f]
   *
   *     where \f$f(x)\f$ and \f$g(y)\f$ are scalar basis functions of \f$U\f$
   *     and \f$V\f$ (*not* any transformations of them, i.e. not their
   *     surface curls, divs etc.), *at least as long as the supports of
   *     \f$f(x)\f$ and \f$g(y)\f$ do not overlap*. (Note that this is
   *     possible---for \f$x \neq y\f$---even for hypersingular operators.)
   *     In this mode, individual blocks of the operator are assembled
   *     differently, depending on whether the bounding box of the supports
   *     of the test basis functions contributing to a given block overlaps
   *     or not with the bounding box of the supports of the trial basis
   *     functions.
   *
   *     Blocks for which these two boxes overlap are assembled as in the
   *     GLOBAL_ASSEMBLY mode, using an arbitrary weak form (possibly one
   *     that does not have the form defined above). For the remaining
   *     blocks, ACA is used to approximate the discretization of the weak
   *     form as defined above with test and trial functions taken as the
   *     restrictions of the the basis functions of \f$U\f$ and \f$V\f$ to
   *     single elements. After the assembly of each such block, linear
   *     combinations of its rows and columns are used to build a low-rank
   *     representation of the block in the original bases of \f$U\f$ and
   *     \f$V\f$.
   *
   *     This mode combines the advantages of LOCAL_ASSEMBLY (fast matrix
   *     construction) and GLOBAL_ASSEMBLY (low memory consumption,
   *     representation in the form of a single H-matrix). However, at
   *     present it cannot be used for Maxwell operators, which require
   *     vector basis functions.
   *
   *  Some integral operators may not support the \p LOCAL_ASSEMBLY and \p
   *  HYBRID_ASSEMBLY modes. How these operators behave when one of these
   *  modes is chosen can be controlled with the #reactionToUnsupportedMode
   *  parameter.
   */
  AcaAssemblyMode mode;

  /** \brief Actions to take when an unsupported assembly mode is detected.
   *
   *  \see reactionToUnsupportedMode. */
  enum ReactionToUnsupportedMode {
    MIN_REACTION,
    IGNORE = MIN_REACTION,
    WARNING,
    ERROR,
    MAX_REACTION = ERROR
  };

  /** \brief Action to take when an unsupported assembly mode is detected.
   *
   *  This parameter controls what happens if an integral operator is
   *  requested to construct its discrete weak form in an ACA mode
   *  (LOCAL_ASSEMBLY or HYBRID_ASSEMBLY) it does not support.
   *
   *  If reactionToUnsupportedMode is set to IGNORE, the operator silently
   *  assembles its discrete weak form in the GLOBAL_ASSEMBLY mode.
   *
   *  It reactionToUnsupportedMode is set to WARNING (default), the operator
   *  emits a warning message to the standard output and assembles its
   *  discrete weak form in the GLOBAL_ASSEMBLY mode.
   *
   *  If reactionToUnsupportedMode is set to ERROR, the operator throws an
   *  exception.
   */
  ReactionToUnsupportedMode reactionToUnsupportedMode;

  /** \brief Recompress ACA matrix after construction?
   *
   *  If true, blocks of H matrices are agglomerated in an attempt to reduce
   *  memory consumption.
   *
   *  Warning: this procedure is not parallelised yet, therefore it may be
   *  slow.
   *
   *  Default value: false. */
  bool recompress;

  /** \brief If true, hierarchical matrix structure will be written in
   *  PostScript format at the end of the assembly procedure.
   *
   *  Default value: false. */
  bool outputPostscript;

  /** \brief Name of the output PostScript file.
   *
   *  \see outputPostscript.
   *
   *  Default value: "aca.ps". */
  std::string outputFname;

  /** \brief Estimate of the magnitude of typical entries of the matrix to be
   *  approximated.
   *
   *  Usually does not need to be changed. Default value: 1. */
  double scaling;

  /** \brief Choice of ACA variant.
   *
   *  If true, the default implementation of the ACA algorithm from the AHMED
   *  library will be used.
   *
   *  If false, an alternative formulation of ACA will be used. This
   *  formulation appears to be more successful in the approximation of
   *  operators whose matrices contain empty blocks (like the double-layer
   *  boundary operator for the Laplace equation) and in general it tries
   *  harder to find a low-rank approximation of a block before falling back
   *  to a dense representation.
   *
   *  Default value: false.
   */
  bool useAhmedAca;

  /** \brief Index of the first block cluster to be approximated using ACA.
   *
   *  This parameter is included to facilitate debugging of ACA algorithms. If
   *  it is set to a nonnegative value, ACA will process block cluster with
   *  index \p firstClusterIndex first.
   *
   *  Default value: -1 (meaning that this parameter is ignored).
   */
  int firstClusterIndex;

  /** \brief Do global assembly before ACA?
   *
   *  \deprecated This parameter is deprecated and should not be used in new
   *code.
   *  If set to true (default), it is ignored. Otherwise it is equivalent to
   *  setting the \p mode parameter to LOCAL_ASSEMBLY.*/
  bool globalAssemblyBeforeCompression;
};

} // namespace Bempp

#endif
