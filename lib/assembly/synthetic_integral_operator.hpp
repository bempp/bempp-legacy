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

#ifndef bempp_synthetic_integral_operator_hpp
#define bempp_synthetic_integral_operator_hpp

#include "../common/common.hpp"

#include "abstract_boundary_operator.hpp"
#include "boundary_operator.hpp"

#include "../common/shared_ptr.hpp"

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForOperators;
/** \endcond */

} // namespace Fiber

namespace Bempp
{

/** \ingroup abstract_boundary_operators
 *  \brief Integral operator with weak form assembled from several pieces.
 *
 *  This class can be used to construct an integral operator assembled in the
 *  following way [1], which in ACA mode may be less time-consuming, though
 *  more memory-consuming. Suppose the formula for the \f$(m, n)\f$th entry of
 *  the weak form of the operator in question can be written as
 *
 *  \f[
 *      A_{mn} = \sum_{\alpha} \int_\Gamma \int_\Sigma
 *      f_{\alpha m}^*(x) \, K_{\alpha}(x, y) \, g_{\alpha n}(y)\,
 *      \mathrm{d}\Gamma(x)\, \mathrm{d}\Sigma(y),
 *  \f]
 *
 *  where \f$f_{\alpha m}(x)\f$ and \f$g_{\alpha n}(y)\f$ are scalar-valued
 *  transformations of test and trial shape functions (e.g. basis functions,
 *  individual components of their surface curls etc.) and \f$K_\alpha(x, y)\f$
 *  kernel functions. Let \f$\mathcal U\f$ and \f$\mathcal V\f$ be two
 *  boundary-element spaces spanned by basis functions \f$u_m : \Gamma \to
 *  \mathbb{R}\f$ and \f$v_n : \Sigma \to \mathbb{R}\f$ such that the support
 *  of each basis function consists of a single element of the respective mesh;
 *  for example, \f$\mathcal U\f$ and \f$\mathcal V\f$ can be the spaces of
 *  piecewise constant or linear, <em>discontinuous</em> functions defined on
 *  \f$\Gamma\f$ and \f$\Sigma\f$. Provided that \f$\mathrm{span}\,f_{\alpha m}
 *  \subset \mathcal U\f$ and \f$\mathrm{span}\,g_{\alpha n} \subset \mathcal
 *  V\f$ for all \f$\alpha\f$, the matrix \f$A\f$ can be rewritten as
 *
 *  \f[
 *      A = \sum_\alpha U_\alpha I_{\mathcal U}^{-1} K_\alpha I_{\mathcal V}^{-1}
 *          V_\alpha,
 *  \f]
 *
 *  with the entries of the individual matrices given by
 *
 *  \f{align*}{
 *      (U_{\alpha})_{mn} &= \int_\Gamma f_{\alpha m}^*(x)\,u_n(x) \,\mathrm{d}\Gamma(x),\\
 *      (I_{\mathcal U})_{mn} &= \int_\Gamma u_m^*(x)\,u_n(x) \,\mathrm{d}\Gamma(x),\\
 *      K_{\alpha mn} &= \int_\Gamma \int_\Sigma u_m^*(x)\,K(x,y)\,v_n(y)
 *          \,\mathrm{d}\Gamma(x)\,\mathrm{d}\Sigma(y),\\
 *      (I_{\mathcal V})_{mn} &= \int_\Sigma v_m^*(y)\,v_n(y) \,\mathrm{d}\Sigma(y),\\
 *      V_{\alpha mn} &= \int_\Sigma v_m^*(y)\,g_{\alpha n}(x) \,\mathrm{d}\Sigma(y).
 *  \f}
 *
 *  All these matrices except \f$K_\alpha\f$ are sparse, and ACA-based assembly
 *  of \f$K_\alpha\f$ can potentially take less time than that of \f$A\f$
 *  because each entry of the former is an integral over a single pair of
 *  elements, while the entries of \f$A\f$ may contain contributions from many
 *  element pairs.
 *
 *  In the current implementation we only handle the case of equal
 *  \f$K_\alpha\f$, i.e. \f$K_\alpha = K\f$ for all \f$\alpha\f$. In addition, we
 *  coalesce \f$U_{\alpha} I_{\mathcal U}^{-1}\f$ and \f$I_{\mathcal V}^{-1}
 *  V_{\alpha}\f$ into single sparse matrices.
 *
 *  [1] This method was proposed by Lars Kielhorn in "On single- and
 *  multi-trace implementations for scattering problems with BETL", 10.
 *  Workshop on Fast Boundary Element Methods in Industrial Applications,
 *  Soellerhaus 27-30 September 2012. */
template <typename BasisFunctionType_, typename ResultType_>
class SyntheticIntegralOperator :
        public AbstractBoundaryOperator<BasisFunctionType_, ResultType_>
{
    typedef AbstractBoundaryOperator<BasisFunctionType_, ResultType_> Base;
public:
    /** \copydoc AbstractBoundaryOperator::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc AbstractBoundaryOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc AbstractBoundaryOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc AbstractBoundaryOperator::QuadratureStrategy */
    typedef typename Base::QuadratureStrategy QuadratureStrategy;
    typedef typename Fiber::LocalAssemblerForOperators<ResultType>
    LocalAssembler;

    /** \brief Constructor.
     *
     *  \param[in] testLocalOps Vector of the operators \f$U_\alpha\f$, as described above.
     *  \param[in] integralOp The operator \f$K\f$, as described above.
     *  \param[in] trialLocalOps Vector of the operators \f$V_\alpha\f$, as described above.
     *  \param[in] label (Optional) operator label.
     *  \param[in] syntheseSymmetry Symmetry flag (see below).
     *
     *  All operators passed as elements of \p testLocalOps and \p trialLocalOps
     *  must be local.
     *
     *  If \p syntheseSymmetry is equal to \p NO_SYMMETRY, either \p testLocalOps or \p
     *  trialLocalOps may be empty. In the former case, \f$A\f$ is assembled as
     *  \f$\sum_\alpha K I_{\mathcal V}^{-1} V_\alpha\f$; in the latter, as
     *  \f$\sum_\alpha U_\alpha I_{\mathcal V}^{-1} K\f$. \p testLocalOps and \p
     *  trialLocalOps must not be empty at the same time.
     *
     *  If \p syntheseSymmetry contains the flag \p SYMMETRIC and/or \p HERMITIAN, \p
     *  trialLocalOps should be left empty and \f$V_\alpha\f$ are taken as the
     *  transposes or Hermitian transposes of \f$U_\alpha\f$.
     */
    SyntheticIntegralOperator(
            const std::vector<BoundaryOperator<BasisFunctionType, ResultType> >& testLocalOps,
            const BoundaryOperator<BasisFunctionType, ResultType>& integralOp,
            const std::vector<BoundaryOperator<BasisFunctionType, ResultType> >& trialLocalOps,
            const std::string& label = "",
            int syntheseSymmetry = NO_SYMMETRY);

    virtual bool isLocal() const;

    /** \brief Get contexts appropriate for construction of internal integral
        operators and auxiliary local operators. */
    static void getContextsForInternalAndAuxiliaryOperators(
            const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
            shared_ptr<const Context<BasisFunctionType, ResultType> >& internalContext,
            shared_ptr<const Context<BasisFunctionType, ResultType> >& auxContext);

protected:
    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormImpl(
            const Context<BasisFunctionType, ResultType>& context) const;

private:
    /** \cond PRIVATE */
    void checkIntegralOperator() const;
    void checkTestLocalOperators() const;
    void checkTrialLocalOperators() const;

private:
    BoundaryOperator<BasisFunctionType, ResultType> m_integralOp;
    std::vector<BoundaryOperator<BasisFunctionType, ResultType> > m_testLocalOps;
    std::vector<BoundaryOperator<BasisFunctionType, ResultType> > m_trialLocalOps;
    int m_syntheseSymmetry;
    /** \endcond */
};

} //namespace Bempp

#endif
