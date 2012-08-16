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

#ifndef bempp_abstract_boundary_operator_hpp
#define bempp_abstract_boundary_operator_hpp

#include "../common/common.hpp"

#include "assembly_options.hpp"
#include "symmetry.hpp"

#include "../common/scalar_traits.hpp"
#include "../common/shared_ptr.hpp"
#include "../fiber/quadrature_strategy.hpp"
#include "../space/space.hpp"

#include <memory>
#include <string>
#include <vector>

namespace arma
{

template <typename eT> class Mat;

}

namespace Bempp
{

class AbstractBoundaryOperatorId;
class Grid;
class GeometryFactory;
template <typename ValueType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class AbstractBoundaryOperatorSum;
template <typename BasisFunctionType, typename ResultType> class Context;
template <typename BasisFunctionType, typename ResultType> class ElementaryAbstractBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename BasisFunctionType, typename ResultType> class ScaledAbstractBoundaryOperator;

/** \ingroup assembly
 *  \brief Abstract boundary operator.
 *
 *  An AbstractBoundaryOperator represents a linear mapping \f$L : X
 *  \to Y\f$ between two function spaces \f$X : S \to K^p\f$
 *  (_domain_) and \f$Y : S \to K^q\f$ (_range_) defined on an
 *  \f$n\f$-dimensional surface \f$S\f$ embedded in an
 *  \f$(n+1)\f$-dimensional domain. \f$K\f$ denotes either the set of
 *  real or complex numbers.
 *
 *  The function assembleWeakForm() can be used to construct the weak
 *  form of the operator.
 *
 *  \tparam BasisFunctionType
 *    Type used to represent components of the functions from the operator's
 *    domain, range and space dual to range.
 *  \tparam ResultType
 *    Type used to represent elements of the weak form of this operator.
 *
 *  Both template parameters can take the following values: \c float, \c
 *  double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>.
 *  Both types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. If \p BasisFunctionType is
 *  set to a complex type, then \p ResultType must be set to the same type.
 */
template <typename BasisFunctionType_, typename ResultType_>
class AbstractBoundaryOperator
{
public:
    /** \brief Type used to represent components of functions on which the operator acts. */
    typedef BasisFunctionType_ BasisFunctionType;
    /** \brief Type used to represent elements of the operator's weak form. */
    typedef ResultType_ ResultType;
    /** \brief Type used to represent coordinates. */
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
    /** \brief Type of the appropriate instantiation of Fiber::QuadratureStrategy. */
    typedef Fiber::QuadratureStrategy<BasisFunctionType, ResultType, GeometryFactory>
    QuadratureStrategy;

    /** @name Construction and destruction
     *  @{ */

    /** \brief Constructor.
     *
     *  \param[in] domain
     *    Function space being the domain of the operator.
     *
     *  \param[in] range
     *    Function space being the range of the operator.
     *
     *  \param[in] dualToRange
     *    Function space dual to the the range of the operator.
     *
     *  \param[in] label
     *    Textual label of the operator (optional, used for debugging).
     *
     *  The spaces \p range and \p dualToRange must be defined on
     *  the same grid.
     */
    AbstractBoundaryOperator(const shared_ptr<const Space<BasisFunctionType> >& domain,
                             const shared_ptr<const Space<BasisFunctionType> >& range,
                             const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
                             const std::string& label = "",
                             const Symmetry symmetry = NO_SYMMETRY);

    /** \brief Destructor. */
    virtual ~AbstractBoundaryOperator();

    /** \brief Identifier.
     *
     *  If the weak form of this operator is cacheable, return a shared pointer
     *  to a valid instance of a subclass of AbstractBoundaryOperatorId that
     *  is guaranteed to be different for all *logically different* abstract
     *  boundary operators.
     *
     *  If the weak form of this operator is not cacheable, return a null shared
     *  pointer. This is the default implementation.
     */
    virtual shared_ptr<const AbstractBoundaryOperatorId> id() const;

    /** @}
     *  @name Spaces
     *  @{ */

    /** \brief Domain.
     *
     *  Return a reference to the function space being the domain of the operator. */
    shared_ptr<const Space<BasisFunctionType> > domain() const;

    /** \brief Range.
     *
     *  Return a reference to the function space being the range of the operator. */
    shared_ptr<const Space<BasisFunctionType> > range() const;

    /** \brief Dual to range.
     *
     *  Return a reference to the function space dual to the range of the operator. */
    shared_ptr<const Space<BasisFunctionType> > dualToRange() const;

    /** @}
     *  @name Label
     *  @{ */

    /** \brief Return the label of the operator. */
    std::string label() const;

    Symmetry symmetry() const;

    /** @}
     *  @name Assembly
     *  @{ */

    /** \brief Return whether this operator is local.
     *
     *  Suppose that an operator \f$A\f$ acting on a function \f$f(x)\f$
     *  produces another function \f$g(x)\f$. We say that \f$A\f$ is local if
     *  the value of \f$g\f$ at any point \f$x\f$ depends only on the values of
     *  \f$f\f$ in an infinitesimal neighbourhood of \f$x\f$.
     *
     *  Multiplicative and differential operators are local and discretization
     *  of their weak forms with finite elements leads to sparse matrices.
     *  Conversely, integral operators are in general non-local and
     *  discretization of their weak forms leads to dense matrices. */
    virtual bool isLocal() const = 0;

    /** \brief Assemble and returns the operator's weak form.
     *
     *  This function constructs a discrete linear operator representing the
     *  matrix \f$L_{jk}\f$ with entries of the form
     *
     *  \f[L_{jk} = \int_S \phi_j L \psi_k,\f]
     *
     *  where \f$L\f$ is the linear operator represented by this object,
     *  \f$S\f$ denotes the surface on which the domain function space \f$X\f$
     *  is defined and which is represented by the grid returned by
     *  <tt>domain.grid()</tt>, \f$\phi_j\f$ is a _trial function_ from the
     *  space \f$Y'\f$ dual to the range of the operator, \f$Y\$, and
     *  \f$\psi_k\f$ is a _test function_ from the domain space \f$X\f$.
     */
    shared_ptr<DiscreteBoundaryOperator<ResultType> > assembleWeakForm(
            const Context<BasisFunctionType, ResultType>& context) const;

protected:
    /** @} */

    /** \brief Given an AssemblyOptions object, construct objects necessary for
     *  subsequent local assembler construction. */
    void collectDataForAssemblerConstruction(
            const AssemblyOptions& options,
            shared_ptr<Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
            shared_ptr<Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
            shared_ptr<GeometryFactory>& testGeometryFactory,
            shared_ptr<GeometryFactory>& trialGeometryFactory,
            shared_ptr<std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
            shared_ptr<std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
            shared_ptr<Fiber::OpenClHandler>& openClHandler,
            bool& cacheSingularIntegrals) const;

    // Probably can be made public
    /** \brief Implementation of the weak-form assembly.
     *
     *  Construct a discrete linear operator representing the matrix \f$L_{jk}\f$
     *  described in assembleWeakForm() and return a shared pointer to it.
     */
    virtual shared_ptr<DiscreteBoundaryOperator<ResultType> >
    assembleWeakFormImpl(const Context<BasisFunctionType, ResultType>& context) const = 0;

private:
    shared_ptr<const Space<BasisFunctionType> > m_domain;
    shared_ptr<const Space<BasisFunctionType> > m_range;
    shared_ptr<const Space<BasisFunctionType> > m_dualToRange;
    std::string m_label;
    Symmetry m_symmetry;
};

} // namespace Bempp

#endif
