%{
#include "grid/grid_view.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/vtk_writer.hpp"
%}

namespace Bempp 
{

%extend GridView 
{
    bool containsEntity(const EntityPointer<0>& ep) const {
        return $self->containsEntity(ep.entity());
    }

    bool containsEntity(const EntityPointer<1>& ep) const {
        return $self->containsEntity(ep.entity());
    }

    bool containsEntity(const EntityPointer<2>& ep) const {
        return $self->containsEntity(ep.entity());
    }

    bool containsEntity(const EntityPointer<3>& ep) const {
        return $self->containsEntity(ep.entity());
    }

    %ignore containsEntity;

    %pythonappend entities %{
        val._parentGrid = self._parentGrid
    %}

    PyObject* entities(int codim) const {
        switch (codim) {
        case 0:
            return SWIG_NewPointerObj(SWIG_as_voidptr($self->entityIterator<0>().release()), 
                                        $descriptor(Bempp::EntityIterator<0>*), SWIG_POINTER_OWN);
        case 1:
            return SWIG_NewPointerObj(SWIG_as_voidptr($self->entityIterator<1>().release()), 
                                        $descriptor(Bempp::EntityIterator<1>*), SWIG_POINTER_OWN);
        case 2:
            return SWIG_NewPointerObj(SWIG_as_voidptr($self->entityIterator<2>().release()), 
                                        $descriptor(Bempp::EntityIterator<2>*), SWIG_POINTER_OWN);
        case 3:
            return SWIG_NewPointerObj(SWIG_as_voidptr($self->entityIterator<3>().release()), 
                                        $descriptor(Bempp::EntityIterator<3>*), SWIG_POINTER_OWN);
        default:
            PyErr_SetString(PyExc_ValueError, "Invalid codimension");
            return NULL;
        }
    }

    // Reference to the parent grid, stored to prevent it from
    // being garbage-collected while this view is alive
    %pythoncode %{
        def parentGrid(self): return self._parentGrid
        parentGrid = property(parentGrid, None, None, "Parent grid")
    %}

    %pythonappend indexSet %{
        val._parentGridView = self
    %}

    %pythonappend vtkWriter %{
        val._parentGridView = self
    %}
}

} // namespace Bempp

%include "grid/grid_view.hpp"
