%{
#include "grid/id_set.hpp"
%}

namespace Bempp 
{

%extend IdSet {
    IdSet::IdType entityId(const EntityPointer<0>& ep) const {
        return $self->entityId(ep.entity());
    }

    IdSet::IdType entityId(const EntityPointer<1>& ep) const {
        return $self->entityId(ep.entity());
    }

    IdSet::IdType entityId(const EntityPointer<2>& ep) const {
        return $self->entityId(ep.entity());
    }

    IdSet::IdType entityId(const EntityPointer<3>& ep) const {
        return $self->entityId(ep.entity());
    }    

    %ignore entityId;

    // Reference to the parent grid, stored to prevent it from
    // being garbage-collected while this view is alive
    %pythoncode %{
        def parentGrid(self): return self._parentGrid
        parentGrid = property(parentGrid, None, None, "Parent grid")
    %}

}

} // namespace Bempp

%include "grid/id_set.hpp"
