#ifndef bempp_aca_approximate_lu_inverse_hpp
#define bempp_aca_approximate_lu_inverse_hpp

#include "discrete_scalar_valued_linear_operator.hpp"

#ifdef WITH_TRILINOS
#include <Thyra_SpmdVectorSpaceBase_decl.hpp>
#endif

#include "ahmed_aux_fwd.hpp"
#include "index_permutation.hpp"

namespace Bempp
{

template <typename ValueType>
class DiscreteAcaScalarValuedLinearOperator;

template <typename ValueType>
class AcaApproximateLuInverse : public DiscreteScalarValuedLinearOperator<ValueType>
{
public:
    /** \brief Construct an approximate LU decomposition of a H-matrix.

    \param[in] fwdOp  Operator represented internally as a H-matrix.
    \param[in] delta  Requested approximation accuracy (M. Bebendorf recommends
                      delta = 0.1). */
    AcaApproximateLuInverse(
            const DiscreteAcaScalarValuedLinearOperator<ValueType>& fwdOp,
            ValueType delta);

    virtual ~AcaApproximateLuInverse();

    virtual void dump() const;

    virtual arma::Mat<ValueType> asMatrix() const;

    virtual unsigned int rowCount() const;
    virtual unsigned int columnCount() const;

    virtual void addBlock(const std::vector<int>& rows,
                          const std::vector<int>& cols,
                          arma::Mat<ValueType>& block) const;

#ifdef WITH_TRILINOS
public:
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > domain() const;
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > range() const;

protected:
    virtual bool opSupportedImpl(Thyra::EOpTransp M_trans) const;
    virtual void applyImpl(
            const Thyra::EOpTransp M_trans,
            const Thyra::MultiVectorBase<ValueType>& X_in,
            const Teuchos::Ptr<Thyra::MultiVectorBase<ValueType> >& Y_inout,
            const ValueType alpha,
            const ValueType beta) const;
#endif

private:
    virtual void applyBuiltInImpl(const TranspositionMode trans,
                                  const arma::Col<ValueType>& x_in,
                                  arma::Col<ValueType>& y_inout,
                                  const ValueType alpha,
                                  const ValueType beta) const;

private:
    typedef bemblcluster<AhmedDofWrapper<ValueType>, AhmedDofWrapper<ValueType> >
    AhmedBemblcluster;

#ifdef WITH_TRILINOS
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType> > m_domainSpace;
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType> > m_rangeSpace;
#else
    unsigned int m_rowCount;
    unsigned int m_columnCount;
#endif

    blcluster* m_blockCluster;
    mblock<ValueType>** m_blocksL;
    mblock<ValueType>** m_blocksU;

    IndexPermutation m_domainPermutation;
    IndexPermutation m_rangePermutation;
};

} // namespace Bempp

#endif
