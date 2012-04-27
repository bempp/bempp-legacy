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
#include "transposition_mode.hpp"

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

namespace Fiber
{

//template <typename BasisValueType, typename ResultType, typename GeometryFactory>
//class LocalAssemblerFactory;

} // namespace Fiber

namespace Bempp
{

class Grid;
class GeometryFactory;
template <typename ValueType> class DiscreteLinearOperator;
template <typename BasisFunctionType, typename ResultType> class ElementaryLinearOperator;
template <typename BasisFunctionType, typename ResultType> class LinearOperatorSuperposition;
template <typename BasisFunctionType, typename ResultType> class GridFunction;

/** \brief "Formal" linear operator.

  This class template is used as a base class for all implementations of
  various types of linear operators, in particular integral operators.

  A LinearOperator represents a linear mapping \f$L : X \to Y\f$, where \f$X :
  S \to K^p\f$ and \f$Y : T \to K^q\f$ are function spaces, with \f$S\f$
  standing for an \f$n\f$-dimensional surface and \f$T\f$ either equal to
  \f$S\f$ or to a \f$(n+1)\f$-dimensional domain in which \f$S\f$ is embedded.
  \f$K\f$ denotes either the set of real or complex numbers.

  The operator is called "formal" because its domain \f$X\f$ is not specified
  yet. The functions assembleWeakForm() and assembleOperator() construct "true"
  linear operators acting on functions from the space passed as the trialSpace
  parameter.

  \tparam ValueType      Type used to represent elements of \f$K\f$. This can be
                         one of: float, double, std::complex<float> and
                         std::complex<double>.
*/
template <typename BasisFunctionType, typename ResultType>
class LinearOperator
{
public:
    typedef Fiber::LocalAssemblerFactory<BasisFunctionType, ResultType, GeometryFactory>
    LocalAssemblerFactory;
    typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;

    LinearOperator(const Space<BasisFunctionType>& testSpace,
                   const Space<BasisFunctionType>& trialSpace);

    LinearOperator(const LinearOperator<BasisFunctionType, ResultType>& other);

    virtual ~LinearOperator();

    /** \brief Assemble the discrete operator */
    void assemble(const LocalAssemblerFactory& factory,
                  const AssemblyOptions& options);

    /** \brief Return \p true if operator is assembled */
    bool isAssembled() const;

    /** \brief Set y_inout := alpha * A * x_in + beta * y_inout, where A is
      this operator. */
    void apply(const TranspositionMode trans,
               const GridFunction<BasisFunctionType, ResultType>& x_in,
               GridFunction<BasisFunctionType, ResultType>& y_inout,
               ResultType alpha, ResultType beta) const;

    /** \brief Return reference to \p DiscreteLinearOperator. */
    const DiscreteLinearOperator<ResultType>& assembledDiscreteLinearOperator() const;

    // Ideas for better names for all methods here are very welcome!!!
    /** \brief Number of components of the functions from the trial space \f$X\f$.

      This is equal to \f$p\f$ in the notation above. */
    virtual int trialComponentCount() const = 0;

    /** \brief Number of components of the functions from the test space \f$Y\f$.

      This is equal to \f$q\f$ in the notation above. */
    virtual int testComponentCount() const = 0;

    /** \brief True if this operator supports the representation \p repr. */
    virtual bool supportsRepresentation(
            AssemblyOptions::Representation repr) const = 0;

    /** \brief Assemble the operator's weak form.

      Construct a discrete linear operator representing the matrix \f$W_{jk}\f$
      whose entries have the form

      \f[W_{jk} = \int_S \phi_j L \psi_k,\f]

      where \f$L\f$ is the linear operator represented by this object, \f$S\f$
      denotes the surface that is the domain of the trial space \f$X\f$ and
      which is represented by the grid returned by trialSpace.grid(),
      \f$\phi_j\f$ is a function from the test space \f$Y\f$ and \f$\psi_k\f$ a
      function from \f$X\f$.

      Note: trialSpace.grid() and testSpace.grid() must return a reference to
      the same object.

      This is the overload intended to be called by the user, who needs to
      provide a local-assembler factory \p factory to be used to construct an
      appropriate local assembler. */
    virtual std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakForm(
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options) const = 0;

    const std::vector<ElementaryLinearOperator<BasisFunctionType, ResultType> const* >&
    localOperators() const;

    const std::vector<ResultType>& multipliers() const;

    const Space<BasisFunctionType>& testSpace() const;
    const Space<BasisFunctionType>& trialSpace() const;

protected:
    void addLocalOperatorsAndMultipliers(
            const std::vector<ElementaryLinearOperator<BasisFunctionType, ResultType> const*>&
            localOperators,
            const std::vector<ResultType>& multipliers);

private:
    std::vector<ElementaryLinearOperator<BasisFunctionType, ResultType> const*>
    m_localOperators;
    std::vector<ResultType> m_multipliers;

    const Space<BasisFunctionType>& m_testSpace;
    const Space<BasisFunctionType>& m_trialSpace;

    std::auto_ptr<const DiscreteLinearOperator<ResultType> > m_discreteOperator;
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
