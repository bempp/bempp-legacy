#ifndef bempp_block_coalescer_hpp
#define bempp_block_coalescer_hpp

#include "../common/common.hpp"
#include "bempp/common/config_ahmed.hpp"

#ifdef WITH_AHMED

#include "aca_options.hpp"
#include "ahmed_aux_fwd.hpp"
#include "../common/boost_scoped_array_fwd.hpp"
#include "../common/boost_shared_array_fwd.hpp"
#include "../common/scalar_traits.hpp"
#include "../common/shared_ptr.hpp"

class Epetra_CrsMatrix;

namespace Bempp
{

template <typename ValueType>
class BlockCoalescer
{
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
    typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
    typedef bbxbemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;
    typedef mblock<typename AhmedTypeTraits<ValueType>::Type> AhmedMblock;

public:
    BlockCoalescer(blcluster* blclustersRoot,
                   blcluster* decomposedBlclustersRoot,
                   const shared_ptr<const Epetra_CrsMatrix>& permutedTestGlobalToFlatLocalMap,
                   const shared_ptr<const Epetra_CrsMatrix>& permutedTrialGlobalToFlatLocalMap,
                   const boost::shared_array<AhmedMblock*>& blocks,
                   const boost::shared_array<AhmedMblock*>& decomposedBlocks,
                   const AcaOptions& acaOptions);

    void coalesceBlock(unsigned index);

private:
    /** \cond PRIVATE */
    void coalesceDenseBlock(unsigned index);
    void coalesceLowRankBlock(unsigned index);

private:
    boost::scoped_array<blcluster*> m_blclusters;
    boost::scoped_array<blcluster*> m_decomposedBlclusters;
    shared_ptr<const Epetra_CrsMatrix> m_permutedTestGlobalToFlatLocalMap;
    shared_ptr<const Epetra_CrsMatrix> m_permutedTrialGlobalToFlatLocalMap;
    boost::shared_array<AhmedMblock*> m_blocks;
    boost::shared_array<AhmedMblock*> m_decomposedBlocks;
    AcaOptions m_acaOptions;
    /** \endcond PRIVATE */
};

} // namespace Bempp

#endif // WITH_AHMED

#endif
