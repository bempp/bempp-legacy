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

#ifndef bempp_discrete_sparse_boundary_operator_hpp
#define bempp_discrete_sparse_boundary_operator_hpp

#include "../common/common.hpp"
#include "bempp/common/config_ahmed.hpp"
#include "bempp/common/config_trilinos.hpp"

#include "discrete_boundary_operator.hpp"

#include "ahmed_aux_fwd.hpp"
#include "symmetry.hpp"
#include "transposition_mode.hpp"

#include "../common/shared_ptr.hpp"
#include "../common/boost_shared_array_fwd.hpp"
#include "../fiber/scalar_traits.hpp"

#ifdef WITH_TRILINOS
#include <Teuchos_RCP.hpp>
#include <Thyra_SpmdVectorSpaceBase_decl.hpp>
/** \cond FORWARD_DECL */
class Epetra_CrsMatrix;
/** \endcond */
#endif

namespace Bempp
{
/** \cond FORWARD_DECL */
class IndexPermutation;
/** \endcond */

/** \ingroup discrete_boundary_operators
 *  \brief Discrete boundary operator stored as a sparse matrix.
 */
template <typename ValueType>
class DiscreteSparseBoundaryOperator :
        public DiscreteBoundaryOperator<ValueType>
{
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
    typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
    typedef bbxbemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;
    typedef mblock<typename AhmedTypeTraits<ValueType>::Type> AhmedMblock;

#ifdef WITH_TRILINOS
public:
    /** \brief Constructor.
     *
     *  \param[in] mat
     *    Sparse matrix that will be represented by the newly
     *    constructed operator. Must not be null.
     *  \param[in] symmetry
     *    Symmetry of the matrix. May be any combination of flags defined
     *    in the Symmetry enumeration type.
     *  \param[in] trans
     *    If different from NO_TRANSPOSE, the discrete operator will represent
     *    a transposed and/or complex-conjugated matrix \p mat. */
    DiscreteSparseBoundaryOperator(
            const shared_ptr<const Epetra_CrsMatrix>& mat,
            int symmetry = NO_SYMMETRY,
            TranspositionMode trans = NO_TRANSPOSE,
            const shared_ptr<AhmedBemBlcluster>& blockCluster =
            shared_ptr<AhmedBemBlcluster>(),
            const shared_ptr<IndexPermutation>& domainPermutation =
            shared_ptr<IndexPermutation>(),
            const shared_ptr<IndexPermutation>& rangePermutation =
            shared_ptr<IndexPermutation>());
#else
    // This class cannot be used without Trilinos
private:
    DiscreteSparseBoundaryOperator();
public:
#endif

    virtual void dump() const;

    virtual arma::Mat<ValueType> asMatrix() const;

    virtual unsigned int rowCount() const;
    virtual unsigned int columnCount() const;

    virtual void addBlock(const std::vector<int>& rows,
                          const std::vector<int>& cols,
                          const ValueType alpha,
                          arma::Mat<ValueType>& block) const;

#ifdef WITH_AHMED
    virtual shared_ptr<const DiscreteBoundaryOperator<ValueType> >
    asDiscreteAcaBoundaryOperator(double eps=-1, int maximumRank=-1, 
        bool interleave=false) const;
#endif

    /** \brief Downcast a shared pointer to a DiscreteBoundaryOperator object to
     *  a shared pointer to a DiscreteSparseBoundaryOperator.
     *
     *  If the object referenced by \p discreteOperator is not in fact a
     *  DiscreteSparseBoundaryOperator, a std::bad_cast exception is thrown. */
    static shared_ptr<const DiscreteSparseBoundaryOperator<ValueType> > castToSparse(
            const shared_ptr<const DiscreteBoundaryOperator<ValueType> >&
            discreteOperator);

#ifdef WITH_TRILINOS
    /** \brief Return a shared pointer to the sparse matrix stored within
     *  this operator.
     *
     *  \note The discrete operator represents the matrix returned by this
     *  function *and possibly transposed and/or complex-conjugated*, depending on
     *  the value returned by transpositionMode(). */
    shared_ptr<const Epetra_CrsMatrix> epetraMatrix() const;
#endif

    /** \brief Return the active sparse matrix transformation.
     *
     *  Indicates whether this operator represents the unmodified sparse matrix
     *  passed in the constructor or its transformation (transposition and/or
     *  conjugation). */
    TranspositionMode transpositionMode() const;

    /** \brief Return the symmetry type of the sparse matrix */
    inline int symmetryMode() const{
        return m_symmetry;
    }

#ifdef WITH_TRILINOS
public:
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > domain() const;
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > range() const;

protected:
    virtual bool opSupportedImpl(Thyra::EOpTransp M_trans) const;
#endif

private:
    /** \cond PRIVATE */
    virtual void applyBuiltInImpl(const TranspositionMode trans,
                                  const arma::Col<ValueType>& x_in,
                                  arma::Col<ValueType>& y_inout,
                                  const ValueType alpha,
                                  const ValueType beta) const;
    bool isTransposed() const;

    // void constructAhmedMatrix(
    //         int* rowOffsets, int* colIndices, double* values,
    //         std::vector<unsigned int>& domain_o2p,
    //         std::vector<unsigned int>& range_p2o,
    //         double eps,
    //         AhmedBemBlcluster* blockCluster,
    //         boost::shared_array<AhmedMblock*>& mblocks,
    //         int& maximumRank) const;
    /** \endcond */

private:
    /** \cond PRIVATE */
#ifdef WITH_TRILINOS
    shared_ptr<const Epetra_CrsMatrix> m_mat;
    int m_symmetry;
    TranspositionMode m_trans;
    shared_ptr<AhmedBemBlcluster> m_blockCluster;
    // o2p
    shared_ptr<IndexPermutation> m_domainPermutation, m_rangePermutation;
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType> > m_domainSpace;
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType> > m_rangeSpace;
#endif
    /** \endcond */
};

} // namespace Bempp

#endif
