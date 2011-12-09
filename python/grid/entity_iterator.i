%{
#include "grid/entity_iterator.hpp"
%}

namespace Bempp
{

%ignore EntityIterator::entity;
%ignore EntityIterator::frozen;

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

} // namespace Bempp

%include "grid/entity_iterator.hpp"

namespace Bempp {
%template(EntityIteratorCodim0) EntityIterator<0>;
%template(EntityIteratorCodim1) EntityIterator<1>;
%template(EntityIteratorCodim2) EntityIterator<2>;
%template(EntityIteratorCodim3) EntityIterator<3>;
} // namespace Bempp