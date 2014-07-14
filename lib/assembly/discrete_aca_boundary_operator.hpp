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
#include "bempp/common/config_ahmed.hpp"

#ifdef WITH_AHMED

#ifndef bempp_discrete_aca_boundary_operator_hpp
#define bempp_discrete_aca_boundary_operator_hpp

#include "../common/common.hpp"

#include "discrete_boundary_operator.hpp"
#include "ahmed_aux_fwd.hpp"
#include "assembly_options.hpp" // actually only ParallelizationOptions are needed
#include "index_permutation.hpp"
#include "symmetry.hpp"
#include "../fiber/scalar_traits.hpp"

#include <iostream>
#include "../common/boost_shared_array_fwd.hpp"

#ifdef WITH_TRILINOS
#include <Teuchos_RCP.hpp>
#include <Thyra_SpmdVectorSpaceBase_decl.hpp>
#endif

namespace Bempp
{

void dumpblcluster(const blcluster* bl, const std::string& indent);

// Forward declarations

/** \cond FORWARD_DECL */
template <typename ValueType> class AcaApproximateLuInverse;
template <typename ValueType> class DiscreteAcaBoundaryOperator;
/** \endcond */

// Global functions

/** \relates DiscreteAcaBoundaryOperator
 *  \brief Add two discrete boundary operators stored as H-matrices.
 *
 *  A std::bad_cast exception is thrown if the input operators can not be
 *  cast to DiscreteAcaBoundaryOperator.
 *
 *  \param[in] op1 First operand.
 *  \param[in] op2 Second operand.
 *  \param[in] eps ??? \todo look into M. Bebendorf's book.
 *  \param[in] maximumRank Maximum rank of blocks that should be considered
 *    low-rank in the H-matrix to be constructed.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the sum of the operands \p op1 and \p op2 stored as a single
 *  H-matrix. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > acaOperatorSum(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2,
        double eps, int maximumRank);

/** \relates DiscreteAcaBoundaryOperator
 *  \brief Multiply the H-matrix representation of a discrete boundary operator
 *  by a scalar and wrap the result in a new discrete boundary operator.
 *
 *  A std::bad_cast exception is thrown if the input operator can not be cast
 *  to DiscreteAcaBoundaryOperator
 *
 *  \param[in] multiplier Scalar multiplier.
 *  \param[in] op Discrete boundary operator to be multiplied.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the operand \p op multiplied by \p multiplier and stored as
 *  a H-matrix. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledAcaOperator(
        const ValueType& multiplier,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op);

/** \relates DiscreteAcaBoundaryOperator
 *  \brief Multiply the H-matrix representation of a discrete boundary operator
 *  by a scalar and wrap the result in a new discrete boundary operator.
 *
 *  A std::bad_cast exception is thrown if the input operator can not be cast
 *  to DiscreteAcaBoundaryOperator
 *
 *  \param[in] op Discrete boundary operator to be multiplied.
 *  \param[in] multiplier Scalar multiplier.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the operand \p op multiplied by \p multiplier and stored as
 *  a H-matrix.  */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledAcaOperator(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        const ValueType& multiplier);

//template <typename ValueType>
//shared_ptr<DiscreteAcaBoundaryOperator<ValueType> > acaOperatorComposition(
//        ValueType multiplier,
//        const shared_ptr<const DiscreteAcaBoundaryOperator<ValueType> >& op1,
//        const shared_ptr<const DiscreteAcaBoundaryOperator<ValueType> >& op2,
//        double eps, int maximumRank);

/** \relates DiscreteAcaBoundaryOperator
 *  \brief LU inverse of a discrete boundary operator stored as a H-matrix.
 *
 *  \param[in] op Discrete boundary operator for which to compute the LU inverse.
 *  \param[in] delta Approximation accuracy of the inverse.
 *
 *  \return A shared pointer to a newly allocated discrete boundary operator
 *  representing the (approximate) LU inverse of \p op and stored as
 *  an (approximate) LU decomposition of \p op. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > acaOperatorApproximateLuInverse(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        double delta);

// class DiscreteAcaBoundaryOperator

/** \ingroup discrete_boundary_operators
 *  \brief Discrete linear operator stored as a H-matrix.
 */
template <typename ValueType>
class DiscreteAcaBoundaryOperator :
        public DiscreteBoundaryOperator<ValueType>
{
    friend class AcaApproximateLuInverse<ValueType>;
    friend shared_ptr<const DiscreteBoundaryOperator<ValueType> > acaOperatorSum<>(
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2,
            double eps, int maximumRank);

    friend shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledAcaOperator<>(
            const ValueType& multiplier,
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op);

    friend shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledAcaOperator<>(
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
            const ValueType& multiplier);

public:
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
    typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
    typedef bbxbemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;
    typedef mblock<typename AhmedTypeTraits<ValueType>::Type> AhmedMblock;
    typedef boost::shared_array<AhmedMblock*> AhmedMblockArray;
    // Unfortunately currently shared_array<T> cannot be converted to
    // shared_array<const T>. So we can't write
    //     typedef boost::shared_array<const AhmedMblock*> AhmedConstMblockArray;
    // and need to use this instead:
    typedef boost::shared_array<AhmedMblock*> AhmedConstMblockArray;

    /** \brief Constructor.
     *
     *  \param[in] rowCount
     *    Number of rows.
     *  \param[in] columnCount
     *    Number of columns.
     *  \param[in] epsUsedInAssembly
     *    The epsilon parameter used during assembly of the H-matrix.
     *  \param[in] maximumRankUsedInAssembly
     *    The limit on block rank used during assembly of the H-matrix.
     *    This is used in particular when generating the approximate H-matrix
     *    LU.
     *  \param[in] symmetry
     *    H-matrix symmetry. Can be any combination of the flags defined in the
     *    Symmetry enumeration type.
     *  \param[in] blockCluster_
     *    Block cluster defining the structure of the H-matrix.
     *  \param[in] blocks_
     *    Array containing the H-matrix blocks.
     *  \param[in] domainPermutation_
     *    Mapping from original to permuted column indices.
     *  \param[in] rangePermutation_
     *    Mapping from original to permuted row indices.
     *  \param[in] parallelizationOptions_
     *    Options determining the maximum number of threads used in
     *    the apply() routine for the H-matrix-vector product.
     *  \param[in] sharedBlocks_
     *    Vector of arrays of mblocks on which this operator implicitly
     *    depends and which therefore must stay alive for the lifetime
     *    of this operator. Useful for constructing ACA operators that
     *    combine mblocks of several other operators.
     *
     *  \note Currently the apply() routine is only parallelized for
     *  non-Hermitian H-matrices.
     */
    DiscreteAcaBoundaryOperator(
            unsigned int rowCount, unsigned int columnCount,
            double epsUsedInAssembly,
            int maximumRankUsedInAssembly,
            int symmetry,
            const shared_ptr<const AhmedBemBlcluster>& blockCluster_,
            const AhmedMblockArray& blocks_,
            const IndexPermutation& domainPermutation_,
            const IndexPermutation& rangePermutation_,
            const ParallelizationOptions& parallelizationOptions_,
            const std::vector<AhmedConstMblockArray>& sharedBlocks_ =
                std::vector<AhmedConstMblockArray>());

    /** \brief Constructor.
     *
     *  \param[in] rowCount
     *    Number of rows.
     *  \param[in] columnCount
     *    Number of columns.
     *  \param[in] maximumRankUsedInAssembly
     *    The limit on block rank used during assembly of the H-matrix.
     *    This is used in particular when generating the approximate H-matrix
     *    LU.
     *  \param[in] symmetry
     *    H-matrix symmetry. Can be any combination of the flags defined in the
     *    Symmetry enumeration type.
     *  \param[in] blockCluster_
     *    Block cluster defining the structure of the H-matrix.
     *  \param[in] blocks_
     *    Array containing the H-matrix blocks.
     *  \param[in] domainPermutation_
     *    Mapping from original to permuted column indices.
     *  \param[in] rangePermutation_
     *    Mapping from original to permuted row indices.
     *  \param[in] parallelizationOptions_
     *    Options determining the maximum number of threads used in
     *    the apply() routine for the H-matrix-vector product.
     *  \param[in] sharedBlocks_
     *    Vector of arrays of mblocks on which this operator implicitly
     *    depends and which therefore must stay alive for the lifetime
     *    of this operator. Useful for constructing ACA operators that
     *    combine mblocks of several other operators.
     *
     *  \note Currently the apply() routine is only parallelized for
     *  non-Hermitian H-matrices.
     *
     *  \deprecated This constructor is deprecated. Use the non-deprecated
     *  constructor. */
    DiscreteAcaBoundaryOperator(
            unsigned int rowCount, unsigned int columnCount,
            int maximumRankUsedInAssembly,
            int symmetry,
            std::unique_ptr<const AhmedBemBlcluster> blockCluster_,
            AhmedMblockArray blocks_,
            const IndexPermutation& domainPermutation_,
            const IndexPermutation& rangePermutation_,
            const ParallelizationOptions& parallelizationOptions_,
            const std::vector<AhmedConstMblockArray>& sharedBlocks_ =
                std::vector<AhmedConstMblockArray>()) BEMPP_DEPRECATED;

    /** \brief Constructor.
     *
     *  \param[in] rowCount
     *    Number of rows.
     *  \param[in] columnCount
     *    Number of columns.
     *  \param[in] maximumRankUsedInAssembly
     *    The limit on block rank used during assembly of the H-matrix.
     *    This is used in particular when generating the approximate H-matrix
     *    LU.
     *  \param[in] symmetry
     *    H-matrix symmetry. Can be any combination of the flags defined in the
     *    Symmetry enumeration type.
     *  \param[in] blockCluster_
     *    Block cluster defining the structure of the H-matrix.
     *  \param[in] blocks_
     *    Array containing the H-matrix blocks.
     *  \param[in] domainPermutation_
     *    Mapping from original to permuted column indices.
     *  \param[in] rangePermutation_
     *    Mapping from original to permuted row indices.
     *  \param[in] parallelizationOptions_
     *    Options determining the maximum number of threads used in
     *    the apply() routine for the H-matrix-vector product.
     *  \param[in] sharedBlocks_
     *    Vector of arrays of mblocks on which this operator implicitly
     *    depends and which therefore must stay alive for the lifetime
     *    of this operator. Useful for constructing ACA operators that
     *    combine mblocks of several other operators.
     *
     *  \note Currently the apply() routine is only parallelized for
     *  non-Hermitian H-matrices.
     *
     *  \deprecated This constructor is deprecated. Use the non-deprecated
     *  constructor.
     */
    BEMPP_DEPRECATED
    DiscreteAcaBoundaryOperator(
            unsigned int rowCount, unsigned int columnCount,
            int maximumRankUsedInAssembly,
            int symmetry,
            const shared_ptr<const AhmedBemBlcluster>& blockCluster_,
            const AhmedMblockArray& blocks_,
            const IndexPermutation& domainPermutation_,
            const IndexPermutation& rangePermutation_,
            const ParallelizationOptions& parallelizationOptions_,
            const std::vector<AhmedConstMblockArray>& sharedBlocks_ =
                std::vector<AhmedConstMblockArray>());

    virtual arma::Mat<ValueType> asMatrix() const;

    virtual unsigned int rowCount() const;
    virtual unsigned int columnCount() const;

    virtual void addBlock(const std::vector<int>& rows,
                          const std::vector<int>& cols,
                          const ValueType alpha,
                          arma::Mat<ValueType>& block) const;

    virtual shared_ptr<const DiscreteBoundaryOperator<ValueType> >
    asDiscreteAcaBoundaryOperator(double eps=-1, int maximumRank=-1,
                                  bool interleave=false) const;

    /** \brief Uncompress all blocks of the H-matrix and store them as dense
     *  matrices.
     *
     *  Sometimes useful for debugging. */
    void makeAllMblocksDense();

    /** \brief Downcast a reference to a DiscreteBoundaryOperator object to
     *  DiscreteAcaBoundaryOperator.
     *
     *  If the object referenced by \p discreteOperator is not in fact a
     *  DiscreteAcaBoundaryOperator, a std::bad_cast exception is thrown. */
    static const DiscreteAcaBoundaryOperator<ValueType>& castToAca(
            const DiscreteBoundaryOperator<ValueType>& discreteOperator);

    /** \brief Downcast a shared pointer to a DiscreteBoundaryOperator object to
     *  a shared pointer to a DiscreteAcaBoundaryOperator.
     *
     *  If the object referenced by \p discreteOperator is not in fact a
     *  DiscreteAcaBoundaryOperator, a std::bad_cast exception is thrown. */
    static shared_ptr<const DiscreteAcaBoundaryOperator<ValueType> > castToAca(
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >&
            discreteOperator);

    /** \brief Return the upper bound for the rank of low-rank mblocks
     *  specified during H-matrix construction. */
    int maximumRank() const;

    /** \brief Return the value of the epsilon parameter specified during
     *  H-matrix construction. */
    double eps() const;

    /** \brief Return the actual maximum rank of low-rank mblocks. */
    int actualMaximumRank() const;

    /** \brief Return a flag describing the symmetry of this operator. */
    int symmetry() const;

    /** \brief Return the block cluster. */
    shared_ptr<const AhmedBemBlcluster> blockCluster() const;

    /** \brief Return the mblocks making up the H-matrix represented by this
     *  operator.
     *
     *  \note This function returns a shared array of pointers to
     *  *non-constant* blocks. However, you must not modify it! This is just
     *  a workaround for AHMED's lack of const-correctness. */
    AhmedMblockArray blocks() const;

    /** \brief Return the number of mblocks making up this operator. */
    size_t blockCount() const;

    /** \brief Return the domain index permutation. */
    const IndexPermutation& domainPermutation() const;

    /** \brief Return the range index permutation. */
    const IndexPermutation& rangePermutation() const;

    /** \brief Return the parallelization options used in the matrix-vector
     *  multiply. */
    const ParallelizationOptions& parallelizationOptions() const;

    /** \brief Return the vector of mblock arrays that this operator implicitly
     *  depends on. */
    std::vector<AhmedConstMblockArray> sharedBlocks() const;

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

    shared_ptr<const AhmedBemBlcluster> m_blockCluster;
    AhmedMblockArray m_blocks;

    IndexPermutation m_domainPermutation;
    IndexPermutation m_rangePermutation;
    ParallelizationOptions m_parallelizationOptions;
    std::vector<AhmedConstMblockArray> m_sharedBlocks;
    /** \endcond */
};

} // namespace Bempp

#endif // WITH_AHMED

#endif
