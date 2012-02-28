#ifndef CONCRETE_ELEMENT_MAPPER_HPP
#define CONCRETE_ELEMENT_MAPPER_HPP

#include "mapper.hpp"

#include <dune/grid/common/mcmgmapper.hh>

namespace Bempp
{

template <typename DuneGridView>
class ConcreteElementMapper : public Mapper<0>
{
private:
    Dune::MultipleCodimMultipleGeomTypeMapper<DuneGridView,
    Dune::MCMGElementLayout > m_dune_mapper_2d;

public:
    ConcreteElementMapper(const DuneGridView& dune_gv) :
        m_dune_mapper_2d(dune_gv,
                         Dune::MCMGElementLayout<DuneGridView::dimension>())
    {}

    /** \brief Map entity to array index. */
    virtual int entityIndex(const Entity<0>& e) const {
        typedef typename DuneGridView::template Codim<0>::Entity DuneEntity;
        typedef ConcreteEntity<0, DuneEntity> ConcEntity;
        const ConcEntity& ce = dynamic_cast<const ConcEntity&>(e);
        return m_dune_mapper_2d.map(ce.duneEntity());
    }

    virtual int subEntityIndex(const Entity<0>& e, int i,
                               unsigned int codimSub) const {
        typedef typename DuneGridView::template Codim<0>::Entity DuneEntity;
        typedef ConcreteEntity<0, DuneEntity> ConcEntity;
        const ConcEntity& ce = dynamic_cast<const ConcEntity&>(e);
        return m_dune_mapper_2d.map(ce.duneEntity(), i, codimSub);
    }
};

}

#endif // CONCRETE_ELEMENT_MAPPER_HPP
