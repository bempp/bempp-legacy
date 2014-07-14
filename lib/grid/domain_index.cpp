#include "domain_index.hpp"

#include "entity_pointer.hpp"
#include "entity.hpp"

#include "../common/acc.hpp"

namespace Bempp
{

int DomainIndex::domain(const Entity<0>& entity) const
{
    std::unique_ptr<EntityPointer<0> > father;
    const Entity<0>* level0Entity = &entity;
    while (level0Entity->hasFather()) {
        father = level0Entity->father();
        level0Entity = &father->entity();
    }
    int index = m_level0IndexSet->entityIndex(*level0Entity);
    return acc(m_domainIndices, index);
}

std::vector<int> DomainIndex::domainIndices() const
{
    return m_domainIndices;
}

} // namespace
