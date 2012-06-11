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

#include "assembly_options.hpp"
#include "symmetry.hpp"
#include "transposition_mode.hpp"

#include "../common/shared_ptr.hpp"
#include "../fiber/local_assembler_factory.hpp"
#include "../space/space.hpp"

#include <boost/mpl/set.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/utility/enable_if.hpp>
#include <memory>
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
template <typename BasisFunctionType, typename ResultType> class LinearOperatorSuperposition;
template <typename BasisFunctionType, typename ResultType> class GridFunction;

/** \ingroup assembly
 *  \brief "Formal" linear operator.
 *
 *  This class template is used as a base class for all implementations of
 *  various types of linear operators, in particular integral operators.
 *
 *  A LinearOperator represents a linear mapping \f$L : X \to Y\f$, where \f$X :
 *  S \to K^p\f$ and \f$Y : T \to K^q\f$ are function spaces, with \f$S\f$
 *  standing for an \f$n\f$-dimensional surface and \f$T\f$ either equal to
 *  \f$S\f$ or to a \f$(n+1)\f$-dimensional domain in which \f$S\f$ is embedded.
 *  \f$K\f$ denotes either the set of real or complex numbers.
 *
 *  The operator is called "formal" because its domain \f$X\f$ is not specified
 *  yet. The functions assembleWeakForm() and assembleDetachedWeakForm()
 *  construct "true" linear operators acting on functions from the space passed
 *  as the trialSpace parameter.
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
     *  \param[in] testSpace
     *    Function space whose elements will be used as test functions during
     *    a subsequent assembly of the weak form of the operator.
     *
     *  \param[in] trialSpace
     *    Function space whose elements will be used as trial functions during
     *    a subsequent assembly of the weak form of the operator.
     *
     *  The objects referenced by \p testSpace and \p trialSpace must continue
     *  to exist at least until the weak form of the operator is assembled.
     */
    LinearOperator(const Space<BasisFunctionType>& testSpace,
                   const Space<BasisFunctionType>& trialSpace);

    /** \brief Copy constructor.
     *
     *  \note If \p other owns an assembled weak form, this weak form is not
     *  copied to the copy-constructed LinearOperator. */
    LinearOperator(const LinearOperator<BasisFunctionType, ResultType>& other);

    /** \brief Destructor. */
    virtual ~LinearOperator();

    /** @}
     *  @name Spaces
     *  @{ */

    /** \brief Test function space.
     *
     *  Return a reference to the function space whose elements are used as
     *  test functions during assembly of the weak form of the operator. */
    const Space<BasisFunctionType>& testSpace() const;

    /** \brief Trial function space.
     *
     *  Return a reference to the function space whose elements are used as
     *  trial functions during assembly of the weak form of the operator. */
    const Space<BasisFunctionType>& trialSpace() const;

    /** \brief Number of components of the functions from the trial space \f$X\f$.

      This is equal to \f$p\f$ in the notation above. */
    virtual int trialComponentCount() const = 0;

    /** \brief Number of components of the functions from the test space \f$Y\f$.

      This is equal to \f$q\f$ in the notation above. */
    virtual int testComponentCount() const = 0;

    /** @}
     *  @name Constituent elementary operators
     *  @{ */

    /** \brief Return a vector of pointers to the elementary linear operators
     *  making up this linear operator.
     *
     *  Each linear operator \f$L\f$ is internally represented by a linear
     *  superposition of <em>elementary</em> linear operators \f$E_i\f$ (see
     *  \ref ElementaryLinearOperator) with weights \f$w_i\f$:
     *
     *  \f[ L = \sum_i w_i E_i. \f]
     *
     *  This function can be used to retrieve the list of the operators \f$E_i\f$.
     *
     *  \sa constituentOperatorWeights.
     */
    virtual std::vector<const ElementaryLinearOperator<BasisFunctionType, ResultType>*>
    constituentOperators() const = 0;

    /** \brief Return weigths of the elementary linear operators
     *  making up this linear operator.
     *
     *  \see constituentOperators.
     */
    virtual std::vector<ResultType> constituentOperatorWeights() const = 0;

    /** @}
     *  @name Assembly
     *  @{ */

    /** \brief Return true if this operator supports the representation \p repr. */
    virtual bool supportsRepresentation(
            AssemblyOptions::Representation repr) const = 0;

    /** \brief Assemble the operator's weak form and store it internally.
     *
     *  This function constructs a discrete linear operator representing the
     *  matrix \f$W_{jk}\f$ with entries of the form
     *
     *  \f[W_{jk} = \int_S \phi_j L \psi_k,\f]
     *
     *  where \f$L\f$ is the linear operator represented by this object, \f$S\f$
     *  denotes the surface that is the domain of the trial space \f$X\f$ and
     *  which is represented by the grid returned by trialSpace.grid(),
     *  \f$\phi_j\f$ is a function from the test space \f$Y\f$ and \f$\psi_k\f$ a
     *  function from \f$X\f$.
     *
     *  The resulting discrete linear operator is stored internally and, if
     *  necessary, can subsequently be accessed via weakForm() or detached
     *  via detachWeakForm(). */
    void assembleWeakForm(const LocalAssemblerFactory& factory,
                          const AssemblyOptions& options,
                          Symmetry symmetry = UNSYMMETRIC);

    /** \brief Assemble and return the operator's weak form.
     *
     * This function constructs the weak form as described in the reference of
     * assembleWeakForm(), but instead of storing the resulting discrete linear
     * operator internally, it transfers its ownership to the caller. */
    std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleDetachedWeakForm(const LocalAssemblerFactory& factory,
                             const AssemblyOptions& options,
                             Symmetry symmetry = UNSYMMETRIC) const;

    /** \brief Return \p true if the operator stores its assembled weak form. */
    bool isWeakFormAssembled() const;

    /** \brief Return a reference to the weak form assembled beforehand.
     *
     * If the weak form has not previously been assembled, a std::runtime_error
     * exception is thrown. */
    const DiscreteLinearOperator<ResultType>& weakForm() const;

    /** \brief Clear the internal pointer to the assembled weak form and
     * transfer its ownership to the caller.
     *
     * \note: Owing to the behaviour of the copy constructor of an auto_ptr, if
     * the value returned by this function is not assigned to anything, the weak
     * form is destroyed and the memory it occupied is freed. */
    std::auto_ptr<DiscreteLinearOperator<ResultType> > detachWeakForm();

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
     *  Construct a discrete linear operator representing the matrix \f$W_{jk}\f$
     *  whose entries have the form
     *
     *  \f[ W_{jk} = \int_S \phi_j L \psi_k, \f]
     *
     *  where \f$L\f$ is the linear operator represented by this object, \f$S\f$
     *  denotes the surface that is the domain of the trial space \f$X\f$ and
     *  which is represented by the grid returned by <tt>trialSpace.grid()</tt>,
     *  \f$\phi_j\f$ is a function from the test space \f$Y\f$ and \f$\psi_k\f$ a
     *  function from \f$X\f$. */
    virtual std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleDetachedWeakFormImpl(const LocalAssemblerFactory& factory,
                                 const AssemblyOptions& options,
                                 Symmetry symmetry) const = 0;

private:
    const Space<BasisFunctionType>& m_testSpace;
    const Space<BasisFunctionType>& m_trialSpace;

    std::auto_ptr<DiscreteLinearOperator<ResultType> > m_weakForm;
};

// Operator overloading

template <typename BasisFunctionType, typename ResultType>
LinearOperatorSuperposition<BasisFunctionType, ResultType> operator+(
        const LinearOperator<BasisFunctionType, ResultType>& op1,
        const LinearOperator<BasisFunctionType, ResultType>& op2);

template <typename BasisFunctionType, typename ResultType>
LinearOperatorSuperposition<BasisFunctionType, ResultType> operator-(
        const LinearOperator<BasisFunctionType, ResultType>& op1,
        const LinearOperator<BasisFunctionType, ResultType>& op2);

// This type machinery is needed to disambiguate between this operator and
// the one taking a LinearOperator and a GridFunction
template <typename BasisFunctionType, typename ResultType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>, std::complex<double> >,
        ScalarType>,
    LinearOperatorSuperposition<BasisFunctionType, ResultType> >::type
operator*(
        const LinearOperator<BasisFunctionType, ResultType>& op,
        const ScalarType& scalar);

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
LinearOperatorSuperposition<BasisFunctionType, ResultType> operator*(
        const ScalarType& scalar,
        const LinearOperator<BasisFunctionType, ResultType>& op);

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
LinearOperatorSuperposition<BasisFunctionType, ResultType> operator/(
        const LinearOperator<BasisFunctionType, ResultType>& op,
        const ScalarType& scalar);

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator*(
        const LinearOperator<BasisFunctionType, ResultType>& op,
        const GridFunction<BasisFunctionType, ResultType>& fun);

} // namespace Bempp

#endif
