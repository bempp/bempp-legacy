%{
#include "grid/grid.hpp"
#include "grid/grid_view.hpp"
#include "grid/entity_iterator.hpp"
%}

namespace Bempp 
{
    
%extend Grid 
{
    bool mark(int refCount, const EntityPointer<0>& ep) {
        return $self->mark(refCount, ep.entity());
    }
    %ignore mark;
    
    int getMark(const Bempp::EntityPointer<0>& ep) const {
        return $self->getMark(ep.entity());
    }
    %ignore getMark;

    %pythonappend leafView %{
        val._parentGrid = self
    %}

    %pythonappend levelView %{
        val._parentGrid = self
    %}

    %pythonappend globalIdSet %{
        val._parentGrid = self
    %}
}

} // namespace Bempp

%include "grid/grid.hpp"
