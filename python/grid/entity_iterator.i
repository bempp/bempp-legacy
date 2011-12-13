// IMPLEMENTATION NOTE:
//
// The object model used to represent entities in Python is different
// than in C++. 
// 
// In C++, we have three template classes -- Entity, EntityPointer and
// EntityIterator -- with EntityIterator derived from
// EntityPointer. Entities are physically stored within entity
// pointers and iterators; the entity() method of EntityPointer and
// Iterator provides a constant reference to the internal entity
// object. Although the reference is constant, the entity is "mutable"
// in the sense that when an iterator is incremented, its internal
// entity is updated (so that any reference stored by the user will
// from then on refer to a different element).
// 
// In Python, such behaviour would surprise the user and would be
// considered as a bug; for example, the following code
//
// # it -> an entity iterator
// all_entities = [e for e in it.entities()]
//
// would return a list of references to the same object, the last
// entity over which the iteration goes. For this reason, the Python
// class Entity manages in fact an instance of the C++ EntityPointer
// class rather than the C++ Entity class. There is no Python class
// named EntityPointer. The Python EntityIterator class exports only
// two methods: __iter__() and next().

%{
#include "grid/entity_iterator.hpp"
%}

%include "entity_iterator_docstrings.i"

namespace Bempp
{

%extend EntityIterator 
{
    %pythoncode %{
        def __iter__(self):
            return self
    %}

    std::auto_ptr<EntityPointer<codim> > next()
    {
        if (!$self->finished()) {
            std::auto_ptr<Bempp::EntityPointer<codim> > result = $self->frozen();
            $self->next();
            return result;
        } 
        else {
            PyErr_SetNone(PyExc_StopIteration);
            return std::auto_ptr<Bempp::EntityPointer<codim> >(NULL);
        }
    }

    %exception next {
        try {
            $action
            if (&*result == NULL) // no exception occurred, but null pointer returned -> stop iteration
            SWIG_fail; // propagate the StopIteration exception
        }
        catch (const std::exception& e) {
            SWIG_exception(SWIG_RuntimeError, e.what());
        }
    }
    
    %ignore next;
}

template<int codim>
class EntityIterator
{
public:
  virtual ~EntityIterator() = 0;
};

%template(EntityIteratorCodim0) EntityIterator<0>;
%template(EntityIteratorCodim1) EntityIterator<1>;
%template(EntityIteratorCodim2) EntityIterator<2>;
%template(EntityIteratorCodim3) EntityIterator<3>;
} // namespace Bempp
