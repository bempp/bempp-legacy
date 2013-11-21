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

#include "bempp/common/config_trilinos.hpp"

#ifndef bempp_discrete_fmm_boundary_operator_hpp
#define bempp_discrete_fmm_boundary_operator_hpp

#include "../common/common.hpp"

#include "../assembly/discrete_boundary_operator.hpp"
#include "../assembly/assembly_options.hpp" // actually only ParallelizationOptions are needed
#include "../assembly/index_permutation.hpp"
#include "../assembly/symmetry.hpp"
#include "../fiber/scalar_traits.hpp"

#include <iostream>
#include "../common/boost_shared_array_fwd.hpp"

#ifdef WITH_TRILINOS
#include <Teuchos_RCP.hpp>
#include <Thyra_SpmdVectorSpaceBase_decl.hpp>
#endif

namespace Bempp
{

// Forward declarations

/** \cond FORWARD_DECL */
template <typename ValueType> class DiscreteFmmBoundaryOperator;
template <typename ValueType> class Octree;
/** \endcond */

// Global functions

/** \relates DiscreteFmmBoundaryOperator
 *  \brief Add two discrete boundary operators stored as H-matrices.
 *
 *  A std::bad_cast exception is thrown if the input operators can not be
 *  cast to DiscreteFmmBoundaryOperator.
 *
 *  \param[in] op1 First operand.
 *  \param[in] op2 Second operand.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the sum of the operands \p op1 and \p op2 stored as a single
 *  H-matrix. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > fmmOperatorSum(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2,
        double eps, int maximumRank);

/** \relates DiscreteFmmBoundaryOperator
 *  \brief Multiply the H-matrix representation of a discrete boundary operator
 *  by a scalar and wrap the result in a new discrete boundary operator.
 *
 *  A std::bad_cast exception is thrown if the input operator can not be cast
 *  to DiscreteFmmBoundaryOperator
 *
 *  \param[in] multiplier Scalar multiplier.
 *  \param[in] op Discrete boundary operator to be multiplied.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the operand \p op multiplied by \p multiplier and stored as
 *  a H-matrix. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledFmmOperator(
        const ValueType& multiplier,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op);

/** \relates DiscreteFmmBoundaryOperator
 *  \brief Multiply the H-matrix representation of a discrete boundary operator
 *  by a scalar and wrap the result in a new discrete boundary operator.
 *
 *  A std::bad_cast exception is thrown if the input operator can not be cast
 *  to DiscreteFmmBoundaryOperator
 *
 *  \param[in] op Discrete boundary operator to be multiplied.
 *  \param[in] multiplier Scalar multiplier.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the operand \p op multiplied by \p multiplier and stored as
 *  a H-matrix.  */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledFmmOperator(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        const ValueType& multiplier);

//template <typename ValueType>
//shared_ptr<DiscreteFmmBoundaryOperator<ValueType> > fmmOperatorComposition(
//        ValueType multiplier,
//        const shared_ptr<const DiscreteFmmBoundaryOperator<ValueType> >& op1,
//        const shared_ptr<const DiscreteFmmBoundaryOperator<ValueType> >& op2,
//        double eps, int maximumRank);

/** \relates DiscreteFmmBoundaryOperator
 *  \brief LU inverse of a discrete boundary operator stored as a H-matrix.
 *
 *  \param[in] op Discrete boundary operator for which to compute the LU inverse.
 *  \param[in] delta Approximation accuracy of the inverse.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the (approximate) LU inverse of \p op and stored as
 *  an (approximate) LU decomposition of \p op. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > fmmOperatorApproximateLuInverse(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        double delta);

// class DiscreteFmmBoundaryOperator

/** \ingroup discrete_boundary_operators
 *  \brief Discrete linear operator stored as a H-matrix.
 */
template <typename ValueType>
class DiscreteFmmBoundaryOperator :
        public DiscreteBoundaryOperator<ValueType>
{
    friend shared_ptr<const DiscreteBoundaryOperator<ValueType> > fmmOperatorSum<>(
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2,
            double eps, int maximumRank);

    friend shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledFmmOperator<>(
            const ValueType& multiplier,
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op);

    friend shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledFmmOperator<>(
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
            const ValueType& multiplier);

public:
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;

    /** \brief Constructor.
     *
     *  \param[in] rowCount
     *    Number of rows.
     *  \param[in] columnCount
     *    Number of columns.
     *  \param[in] symmetry
     *    H-matrix symmetry. Can be any combination of the flags defined in the
     *    Symmetry enumeration type.
     *
     */
    DiscreteFmmBoundaryOperator(
			unsigned int rowCount, unsigned int columnCount,
			const shared_ptr<Octree<ValueType> > &octree,
			int symmetry);

    virtual arma::Mat<ValueType> asMatrix() const;

    virtual unsigned int rowCount() const;
    virtual unsigned int columnCount() const;

    virtual void addBlock(const std::vector<int>& rows,
                          const std::vector<int>& cols,
                          const ValueType alpha,
                          arma::Mat<ValueType>& block) const;

    virtual shared_ptr<const DiscreteBoundaryOperator<ValueType> >
    asDiscreteFmmBoundaryOperator(double eps=-1, int maximumRank=-1,
                                  bool interleave=false) const;

    /** \brief Downcast a reference to a DiscreteBoundaryOperator object to
     *  DiscreteFmmBoundaryOperator.
     *
     *  If the object referenced by \p discreteOperator is not in fact a
     *  DiscreteFmmBoundaryOperator, a std::bad_cast exception is thrown. */
    static const DiscreteFmmBoundaryOperator<ValueType>& castToFmm(
            const DiscreteBoundaryOperator<ValueType>& discreteOperator);

    /** \brief Downcast a shared pointer to a DiscreteBoundaryOperator object to
     *  a shared pointer to a DiscreteFmmBoundaryOperator.
     *
     *  If the object referenced by \p discreteOperator is not in fact a
     *  DiscreteFmmBoundaryOperator, a std::bad_cast exception is thrown. */
    static shared_ptr<const DiscreteFmmBoundaryOperator<ValueType> > castToFmm(
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >&
            discreteOperator);

    /** \brief Return the upper bound for the rank of low-rank mblocks
     *  specified during H-matrix construction. */
    int maximumRank() const;

    /** \brief Return the value of the epsilon parameter specified during
     *  H-matrix construction. */
    double eps() const;

    /** \brief Return a flag describing the symmetry of this operator. */
    int symmetry() const;

    /** \brief Return the domain index permutation. */
    const IndexPermutation& domainPermutation() const;

    /** \brief Return the range index permutation. */
    const IndexPermutation& rangePermutation() const;

    /** \brief Return the parallelization options used in the matrix-vector
     *  multiply. */
    const ParallelizationOptions& parallelizationOptions() const;

#ifdef WITH_TRILINOS
public:
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > domain() const;
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > range() const;

protected:
    virtual bool opSupportedImpl(Thyra::EOpTransp M_trans) const;
#endif

private:
    virtual void applyBuiltInImpl(const TranspositionMode trans,
                                  const arma::Col<ValueType>& x_in,
                                  arma::Col<ValueType>& y_inout,
                                  const ValueType alpha,
                                  const ValueType beta) const;

private:
    /** \cond PRIVATE */
#ifdef WITH_TRILINOS
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType> > m_domainSpace;
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType> > m_rangeSpace;
#else
    unsigned int m_rowCount;
    unsigned int m_columnCount;
#endif
    double m_eps;
    int m_maximumRank; // used by the approximate-LU preconditioner
    int m_symmetry;

    //IndexPermutation m_domainPermutation;
    //IndexPermutation m_rangePermutation;
    ParallelizationOptions m_parallelizationOptions;
    shared_ptr<Octree<ValueType> > m_octree;
    /** \endcond */
};

} // namespace Bempp

#endif
