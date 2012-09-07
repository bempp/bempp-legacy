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

/** \cond FORWARD_DECL */
template <typename eT> class Mat;
/** \endcond */

}

namespace Bempp
{

/** \cond FORWARD_DECL */
class AbstractBoundaryOperatorId;
class Grid;
class GeometryFactory;
template <typename ValueType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class AbstractBoundaryOperatorSum;
template <typename BasisFunctionType, typename ResultType> class Context;
template <typename BasisFunctionType, typename ResultType> class ElementaryAbstractBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename BasisFunctionType, typename ResultType> class ScaledAbstractBoundaryOperator;
/** \endcond */

/** \ingroup abstract_boundary_operators
 *  \brief Abstract (non-discretized) boundary operator.
 *
 *  An AbstractBoundaryOperator represents a linear mapping \f$L : X
 *  \to Y\f$ between two function spaces \f$X : S \to K^p\f$
 *  (_domain_) and \f$Y : T \to Q^q\f$ (_range_) defined on
 *  \f$n\f$-dimensional surfaces \f$S\f$ and \f$T\f$ embedded in an
 *  \f$(n+1)\f$-dimensional domain. Each of the symbols \f$K\f$ and
 *  \f$Q\f$ can stand either for the set of real or complex numbers.
 *  The surfaces \f$S\f$ and \f$T\f$ may be equal.
 *
 *  The function assembleWeakForm() can be used to construct the weak
 *  form of the operator.
 *
 *  \tparam BasisFunctionType_
 *    Type of the values of the (components of the) basis functions into
 *    which functions acted upon by the operator are expanded.
 *  \tparam ResultType_
 *    Type used to represent elements of the weak form of the operator.
 *
 *  Both template parameters can take the following values: \c float, \c
 *  double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>.
 *  Both types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. If \p BasisFunctionType_ is
 *  set to a complex type, then \p ResultType_ must be set to the same type.
 */
template <typename BasisFunctionType_, typename ResultType_>
class AbstractBoundaryOperator
{
public:
    /** \brief Type of the values of the (components of the) basis functions into
     *  which functions acted upon by the operator are expanded. */
    typedef BasisFunctionType_ BasisFunctionType;
    /** \brief Type used to represent elements of the weak form of the operator. */
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
     *  \param[in] range
     *    Function space being the range of the operator.
     *  \param[in] dualToRange
     *    Function space dual to the the range of the operator.
     *  \param[in] label
     *    Textual label of the operator. If empty, a unique label is generated
     *    automatically.
     *  \param[in] symmetry
     *    Symmetry of the weak form of the operator. Can be any combination of
     *    the flags defined in the enumeration type Symmetry.
     *
     *  None of the shared pointers may be null and the spaces \p range and \p
     *  dualToRange must be defined on the same grid, otherwise an exception is
     *  thrown. */
    AbstractBoundaryOperator(
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
            const std::string& label,
            int symmetry);

    /** \brief Destructor. */
    virtual ~AbstractBoundaryOperator();

    /** \brief Return the identifier of this operator.
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
     *  Return a shared pointer to the function space being the domain of
     *  this operator. */
    shared_ptr<const Space<BasisFunctionType> > domain() const;

    /** \brief Range.
     *
     *  Return a shared pointer to the function space being the range of
     *  this operator. */
    shared_ptr<const Space<BasisFunctionType> > range() const;

    /** \brief Dual to range.
     *
     *  Return a shared pointer to the function space dual to the range of
     *  this operator. */
    shared_ptr<const Space<BasisFunctionType> > dualToRange() const;

    /** @}
     *  @name Other attributes
     *  @{ */

    /** \brief Return the label of the operator. */
    std::string label() const;

    /** \brief Return the symmetry properties of the operator.
     *
     *  The returned value should be treated as a bitwise combination
     *  of the flags defined in the Symmetry enumeration type. */
    int symmetry() const;

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

    /** @}
     *  @name Assembly
     *  @{ */

    /** \brief Assemble and return the operator's weak form.
     *
     *  This function constructs a discrete linear operator representing the
     *  matrix \f$L_{jk}\f$ with entries of the form
     *
     *  \f[L_{jk} = \int_S \phi_j L \psi_k,\f]
     *
     *  where \f$L\f$ is the linear operator represented by this object,
     *  \f$S\f$ denotes the surface on which the domain function space \f$X\f$
     *  is defined and which is represented by the grid returned by
     *  <tt>domain.grid()</tt>, \f$\phi_j\f$ is a _test function_ from the
     *  space \f$Y'\f$ dual to the range of the operator, \f$Y\f$, and
     *  \f$\psi_k\f$ is a _trial function_ from the domain space \f$X\f$.
     */
    shared_ptr<DiscreteBoundaryOperator<ResultType_> > assembleWeakForm(
            const Context<BasisFunctionType_, ResultType_>& context) const;

protected:

    /** \brief Given an AssemblyOptions object, construct objects necessary for
     *  subsequent local assembler construction. */
    void collectDataForAssemblerConstruction(
            const AssemblyOptions& options,
            shared_ptr<Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
            shared_ptr<Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
            shared_ptr<GeometryFactory>& testGeometryFactory,
            shared_ptr<GeometryFactory>& trialGeometryFactory,
            shared_ptr<std::vector<const Fiber::Basis<BasisFunctionType_>*> >& testBases,
            shared_ptr<std::vector<const Fiber::Basis<BasisFunctionType_>*> >& trialBases,
            shared_ptr<Fiber::OpenClHandler>& openClHandler,
            bool& cacheSingularIntegrals) const;

    /** \brief Assemble and return the operator's weak form.
     *
     *  This virtual function is invoked by assembleWeakForm() to do the actual
     *  work. */
    virtual shared_ptr<DiscreteBoundaryOperator<ResultType> >
    assembleWeakFormImpl(const Context<BasisFunctionType, ResultType>& context) const = 0;

    /** @} */

private:
    /** \cond PRIVATE */
    shared_ptr<const Space<BasisFunctionType> > m_domain;
    shared_ptr<const Space<BasisFunctionType> > m_range;
    shared_ptr<const Space<BasisFunctionType> > m_dualToRange;
    std::string m_label;
    int m_symmetry;
    /** \endcond */
};

} // namespace Bempp

#endif
