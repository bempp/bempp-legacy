%{
#include "grid/index_set.hpp"
%}

%include "index_set_docstrings.i"

namespace Bempp
{

%extend IndexSet {
    IndexSet::IndexType entityIndex(const EntityPointer<0>& ep) const {
        return $self->entityIndex(ep.entity());
    }

    IndexSet::IndexType entityIndex(const EntityPointer<1>& ep) const {
        return $self->entityIndex(ep.entity());
    }

    IndexSet::IndexType entityIndex(const EntityPointer<2>& ep) const {
        return $self->entityIndex(ep.entity());
    }

    IndexSet::IndexType entityIndex(const EntityPointer<3>& ep) const {
        return $self->entityIndex(ep.entity());
    }    

    IndexSet::IndexType subEntityIndex(const EntityPointer<0>& ep, int i, unsigned int codimSub) const {
        return $self->subEntityIndex(ep.entity(), i, codimSub);
    }    

    %ignore entityIndex;
    %ignore subEntityIndex;

    // Reference to the parent grid view, stored to prevent it from
    // being garbage-collected while this index set is alive
    %pythoncode %{
        def parentGridView(self): return self._parentGridView
        parentGridView = property(parentGridView, None, None, "Parent grid view")
    %}
}

} // namespace Bempp

%include "grid/index_set.hpp"
