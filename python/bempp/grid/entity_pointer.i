%{
#include "grid/entity_pointer.hpp"
%}

%include "entity_pointer_docstrings.i"

namespace Bempp 
{

%ignore EntityPointer::entity;

template<int codim> class EntityIterator;

%extend EntityPointer 
{
    /* This might be the simplest solution -- to overload EntityPointer's operator->().
    Then it wouldn't be necessary to forward by hand all the method calls below. 
    */
    /*
    const Entity<codim>* operator->() const
    {
    return &$self->entity();
    }*/

    int level() const {
        return $self->entity().level();
    }

    %pythonappend geometry %{
        val._parentEntity = self
    %}

    const Geometry& geometry() const {
        return $self->entity().geometry();
    }
    
    GeometryType type() const {
        return $self->entity().type();
    }

    // Reference to the parent grid, stored to prevent it from
    // being garbage-collected while this entity pointer is alive
    %pythoncode %{
        def parentGrid(self): return self._parentGrid
        parentGrid = property(parentGrid, None, None, "Parent grid")
    %}        
}

%extend EntityPointer<0> 
{
    int subEntityCount(int codim) const {
        switch (codim) {
        case 1:
            return $self->entity().subEntityCount<1>();
        case 2:
            return $self->entity().subEntityCount<2>();
        case 3:
            return $self->entity().subEntityCount<3>();
        default:
            // We assume that the behaviour for all invalid dimensions is the same.
            return $self->entity().subEntityCount<-1>();
        }
    }

    %pythonappend subEntities %{
        val._parentEntity = self # there is no need to expose this attribute to the user
        val._parentGrid = self._parentGrid
    %}
    
    PyObject* subEntities(int codim) const {
        switch (codim)
        {
        case 1:
            return SWIG_NewPointerObj(SWIG_as_voidptr($self->entity().subEntityIterator<1>().release()), 
                                        $descriptor(Bempp::EntityIterator<1>*), SWIG_POINTER_OWN);
        case 2:
            return SWIG_NewPointerObj(SWIG_as_voidptr($self->entity().subEntityIterator<2>().release()), 
                                        $descriptor(Bempp::EntityIterator<2>*), SWIG_POINTER_OWN);
        case 3:
            return SWIG_NewPointerObj(SWIG_as_voidptr($self->entity().subEntityIterator<3>().release()), 
                                        $descriptor(Bempp::EntityIterator<3>*), SWIG_POINTER_OWN);
        default:
            PyErr_SetString(PyExc_ValueError, "Invalid codimension");
            return NULL;
        }
    }
    
    %pythonappend father %{
        if val is not None:
            val._parentGrid = self._parentGrid
    %}
    
    std::auto_ptr<EntityPointer<0> > father() const {
        return $self->entity().father();
    }
    
    bool hasFather() const {
        return $self->entity().hasFather();
    }
    
    bool isLeaf() const {
        return $self->entity().isLeaf();
    }
    
    bool isRegular() const {
        return $self->entity().isRegular();
    }
    
    %pythonappend sons %{
        val._parentEntity = self # there is no need to expose this attribute to the user
        val._parentGrid = self._parentGrid
    %}
    
    std::auto_ptr<EntityIterator<0> > sons(int maxLevel) const {
        return $self->entity().sonIterator(maxLevel);
    }
    
    bool isNew() const {
        return $self->entity().isNew();
    }
    
    bool mightVanish() const {
        return $self->entity().mightVanish();
    }
}

} // namespace Bempp

%include "grid/entity_pointer.hpp"

namespace Bempp 
{
// Note that in Python we rename the C++ EntityPointer to Entity!
%template(EntityCodim0) EntityPointer<0>;
%template(EntityCodim1) EntityPointer<1>;
%template(EntityCodim2) EntityPointer<2>;
%template(EntityCodim3) EntityPointer<3>;
} // namespace Bempp
