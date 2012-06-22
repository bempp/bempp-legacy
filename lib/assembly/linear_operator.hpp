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

#ifndef bempp_linear_operator_hpp
#define bempp_linear_operator_hpp

#include "../common/common.hpp"

#include "assembly_options.hpp"
#include "symmetry.hpp"
#include "transposition_mode.hpp"

#include "../common/scalar_traits.hpp"
#include "../common/shared_ptr.hpp"
#include "../fiber/local_assembler_factory.hpp"
#include "../space/space.hpp"

#include <boost/mpl/set.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/utility/enable_if.hpp>
#include <memory>
#include <string>
#include <vector>

namespace arma
{

template <typename eT> class Mat;

}

namespace Bempp
{

class Grid;
class GeometryFactory;
template <typename ValueType> class DiscreteLinearOperator;
template <typename BasisFunctionType, typename ResultType> class ElementaryLinearOperator;
template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename BasisFunctionType_, typename ResultType_> class ScaledLinearOperator;
template <typename BasisFunctionType_, typename ResultType_> class LinearOperatorSum;

/** \ingroup assembly
 *  \brief Boundary linear operator.
 *
 *  A BoundaryOperator represents a linear mapping \f$L : X \to Y\f$ between
 *  two function spaces \f$X : S \to K^p\f$ (_domain_) and \f$Y : S \to K^q\f$ (_range_) defined on
 *  an \f$n\f$-dimensional surface \f$S\f$ embedded in an
 *  \f$(n+1)\f$-dimensional domain. \f$K\f$ denotes either the set of real or
 *  complex numbers.
 *
 *  The functions assembleWeakForm() and assembleDetachedWeakForm() can be used
 *  to calculate the weak form of the operator. The weak form constructed with
 *  the former function is stored internally in the BoundaryOperator object,
 *  which can subsequently be used as a typical linear operator (i.e. act on
 *  functions defined on the surface \f$S\f$, represented as GridFunction
 *  objects). The weak form constructed with the latter function is not stored
 *  in BoundaryOperator, but its ownership is passed directly to the caller.
 *
 *  \tparam BasisFunctionType
 *    Type used to represent components of the test and trial functions.
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
class LinearOperator
{
public:
    /** \brief Type used to represent components of the test and trial functions. */
    typedef BasisFunctionType_ BasisFunctionType;
    /** \brief Type used to represent elements of the operator's weak form. */
    typedef ResultType_ ResultType;
    /** \brief Type used to represent coordinates. */
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
    /** \brief Type of the appropriate instantiation of Fiber::LocalAssemblerFactory. */
    typedef Fiber::LocalAssemblerFactory<BasisFunctionType, ResultType, GeometryFactory>
    LocalAssemblerFactory;

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
     *  The objects referenced by \p domain, \p range and \p dualToRange must
     *  continue to exist at least until the weak form of the operator is
     *  assembled. The spaces \p range and \p dualToRange must be defined on
     *  the same grid.
     */
    LinearOperator(const Space<BasisFunctionType>& domain,
                   const Space<BasisFunctionType>& range,
                   const Space<BasisFunctionType>& dualToRange,
                   const std::string& label = "");

    // Default "shallow" copy constructor is used (thus the internal pointer
    // to the weak form is shared among all copies). To make a deep copy of
    // a LinearOperator, use the deepCopy() method (not yet written).
    // LinearOperator(const LinearOperator<BasisFunctionType, ResultType>& other);

    /** \brief Destructor. */
    virtual ~LinearOperator();

    /** \brief Clone operator.
     *
     *  Returns an auto pointer to a shallow copy (i.e. sharing a
     *  reference-counted pointer to the weak form) of this operator
     *  allocated on the heap.
     */
    virtual std::auto_ptr<LinearOperator> clone() const = 0;

    /** @}
     *  @name Spaces
     *  @{ */

    /** \brief Domain.
     *
     *  Return a reference to the function space being the domain of the operator. */
    const Space<BasisFunctionType>& domain() const;

    /** \brief Range.
     *
     *  Return a reference to the function space being the range of the operator. */
    const Space<BasisFunctionType>& range() const;

    /** \brief Dual to range.
     *
     *  Return a reference to the function space dual to the range of the operator. */
    const Space<BasisFunctionType>& dualToRange() const;

    /** @}
     *  @name Label
     *  @{ */

    /** \brief Return the label of the operator. */
    std::string label() const;

    /** \brief Set the label of the operator. */
    void setLabel(const std::string& newLabel);

    /** @}
     *  @name Assembly
     *  @{ */

    /** \brief Return true if this operator supports the representation \p repr. */
    virtual bool supportsRepresentation(
            AssemblyOptions::Representation repr) const = 0;

    //*  \param[in] force
    //*    If true (default), the weak form will be reassembled even if an older one already exists. If false, any existing a weak form
    /** \brief Assemble the operator's weak form and store it internally.
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
     *
     *  The resulting discrete linear operator is stored internally. It can
     *  subsequently be accessed via weakForm() or, if necessary, reset via
     *  resetWeakForm().
     *
     *  Note that the weak form is assembled anew even if the function
     *  has already been called before. In order to avoid a repeated assembly,
     *  simply write
     *
     *  \code{cpp}
     *  if (!op.isWeakFormAssembled())
     *      op.assembleWeakForm(...);
     *  \endcode
     */
    void assembleWeakForm(const LocalAssemblerFactory& factory,
                          const AssemblyOptions& options,
                          Symmetry symmetry = NO_SYMMETRY/*,
                          bool force = true*/);

    /** \brief Return \p true if the operator stores its assembled weak form. */
    bool isWeakFormAssembled() const;

    /** \brief Return a shared pointer to the weak form assembled beforehand.
     *
     * If the weak form has not previously been assembled, a std::runtime_error
     * exception is thrown. */
    shared_ptr<const DiscreteLinearOperator<ResultType> > weakForm() const;

    /** \brief Reset the internal shared pointer to the weak form.
     *
     *  Note that other objects (e.g. composite operators) may also store shared
     *  pointers to this weak form. Its memory will only be freed after the
     *  last such shared pointer is reset. */
    void resetWeakForm();

    /** @}
     *  @name Action
     *  @{ */

    /** \brief Set <tt>y_inout := alpha * A * x_in + beta * y_inout</tt>, where
     *  \c A is this operator. */
    void apply(const TranspositionMode trans,
               const GridFunction<BasisFunctionType, ResultType>& x_in,
               GridFunction<BasisFunctionType, ResultType>& y_inout,
               ResultType alpha, ResultType beta) const;

    /** @} */

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

    /** \brief Implementation of the weak-form assembly.
     *
     *  Construct a discrete linear operator representing the matrix \f$L_{jk}\f$
     *  described in assembleWeakForm() and return a shared pointer to it.
     */
    virtual shared_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakFormImpl(const LocalAssemblerFactory& factory,
                         const AssemblyOptions& options,
                         Symmetry symmetry) = 0;

private:
    const Space<BasisFunctionType>& m_domain;
    const Space<BasisFunctionType>& m_range;
    const Space<BasisFunctionType>& m_dualToRange;
    std::string m_label;

    shared_ptr<DiscreteLinearOperator<ResultType> > m_weakForm;
};

} // namespace Bempp

// The includes below are not strictly needed, but convenient for users
// employing the arithmetic operator overloads
#include "scaled_linear_operator.hpp"
#include "linear_operator_sum.hpp"
// #include "linear_operator_composition.hpp"

// Operator overloading

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
LinearOperatorSum<BasisFunctionType, ResultType> operator+(
        const LinearOperator<BasisFunctionType, ResultType>& op1,
        const LinearOperator<BasisFunctionType, ResultType>& op2);

template <typename BasisFunctionType, typename ResultType>
LinearOperatorSum<BasisFunctionType, ResultType> operator-(
        const LinearOperator<BasisFunctionType, ResultType>& op1,
        const LinearOperator<BasisFunctionType, ResultType>& op2);

// This type machinery is needed to disambiguate between this operator and
// the one taking a LinearOperator and a GridFunction
template <typename BasisFunctionType, typename ResultType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>, std::complex<double> >,
        ScalarType>,
    ScaledLinearOperator<BasisFunctionType, ResultType> >::type
operator*(
        const LinearOperator<BasisFunctionType, ResultType>& op,
        const ScalarType& scalar);

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
ScaledLinearOperator<BasisFunctionType, ResultType> operator*(
        const ScalarType& scalar,
        const LinearOperator<BasisFunctionType, ResultType>& op);

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
ScaledLinearOperator<BasisFunctionType, ResultType> operator/(
        const LinearOperator<BasisFunctionType, ResultType>& op,
        const ScalarType& scalar);

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator*(
        const LinearOperator<BasisFunctionType, ResultType>& op,
        const GridFunction<BasisFunctionType, ResultType>& fun);

//template <typename BasisFunctionType, typename ResultType>
//LinearOperatorComposition<BasisFunctionType, ResultType> operator*(
//        const LinearOperator<BasisFunctionType, ResultType>& op1,
//        const LinearOperator<BasisFunctionType, ResultType>& op2);

} // namespace Bempp

#endif
