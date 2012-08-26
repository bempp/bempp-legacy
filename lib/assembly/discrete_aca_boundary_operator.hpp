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

#ifndef bempp_discrete_aca_boundary_operator_hpp
#define bempp_discrete_aca_boundary_operator_hpp

#include "../common/common.hpp"

#include "discrete_boundary_operator.hpp"
#include "ahmed_aux_fwd.hpp"
#include "index_permutation.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../common/not_implemented_error.hpp"

#include <iostream>
#include "../common/boost_shared_array_fwd.hpp"

#ifdef WITH_TRILINOS
#include <Teuchos_RCP.hpp>
#include <Thyra_SpmdVectorSpaceBase_decl.hpp>
#endif

namespace Bempp {

template <typename ValueType> class AcaApproximateLuInverse;

/** \ingroup discrete_boundary_operators
 *  \brief Discrete linear operator stored as a hierarchical (compressed) matrix.
 */
template <typename ValueType>
class DiscreteAcaBoundaryOperator :
        public DiscreteBoundaryOperator<ValueType>
{
    friend class AcaApproximateLuInverse<ValueType>;

public:
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
    typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
    typedef bemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;
    typedef mblock<typename AhmedTypeTraits<ValueType>::Type> AhmedMblock;

    DiscreteAcaBoundaryOperator(
            unsigned int rowCount, unsigned int columnCount,
            int maximumRank,
            bool symmetric,
            std::auto_ptr<AhmedBemBlcluster> blockCluster,
            boost::shared_array<AhmedMblock*> blocks,
            const IndexPermutation& domainPermutation,
            const IndexPermutation& rangePermutation);

    virtual arma::Mat<ValueType> asMatrix() const;

    virtual unsigned int rowCount() const;
    virtual unsigned int columnCount() const;

    virtual void addBlock(const std::vector<int>& rows,
                          const std::vector<int>& cols,
                          const ValueType alpha,
                          arma::Mat<ValueType>& block) const;

    void makeAllMblocksDense(); // for debugging

    static const DiscreteAcaBoundaryOperator<ValueType>& castToAca(
            const DiscreteBoundaryOperator<ValueType>& discreteOperator);

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
#ifdef WITH_TRILINOS
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType> > m_domainSpace;
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType> > m_rangeSpace;
#else
    unsigned int m_rowCount;
    unsigned int m_columnCount;
#endif
    int m_maximumRank; // used by the approximate-LU preconditioner
    bool m_symmetric;

    std::auto_ptr<AhmedBemBlcluster> m_blockCluster;
    boost::shared_array<AhmedMblock*> m_blocks;

    IndexPermutation m_domainPermutation;
    IndexPermutation m_rangePermutation;
};


} // namespace Bempp

#endif