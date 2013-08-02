%{
#include "grid/grid_view.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/vtk_writer.hpp"
%}

%include "grid_view_docstrings.i"

// Handle the enum Dune::VTK::OutputType like a string
%typemap(in) Dune::VTK::DataMode
{
    if (!PyString_Check($input))
    {
        PyErr_SetString(PyExc_TypeError, "in method '$symname', argument $argnum: expected a string");
        SWIG_fail;
    }
    const std::string s(PyString_AsString($input));
    if (s == "conforming")
        $1 = Dune::VTK::conforming;
    else if (s == "nonconforming")
        $1 = Dune::VTK::nonconforming;
    else
    {
        PyErr_SetString(PyExc_ValueError, "in method '$symname', argument $argnum: expected one of 'conforming' or 'nonconforming'");
        SWIG_fail;
    }
}

%apply arma::Mat<double>& ARGOUT_MAT {
       arma::Mat<double>& vertices
};

%apply arma::Mat<int>& ARGOUT_MAT {
  arma::Mat<int>& elementCorners
};

%apply arma::Mat<char>& ARGOUT_MAT {
  arma::Mat<char>& auxData
};

%apply arma::Mat<char>& ARGOUT_MAT {
  arma::Mat<char>& auxData
};

%apply std::vector<int>& ARGOUT_VEC {
    std::vector<int>& domainIndices
};

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

    %ignore elementMapper;
    %ignore reverseElementMapper;

    %rename(_getRawElementDataNoElementIndices) getRawElementData(
        arma::Mat<double>&, arma::Mat<int>&, arma::Mat<char>&) const;
    %rename(_getRawElementDataWithElementIndices) getRawElementData(
        arma::Mat<double>&, arma::Mat<int>&, arma::Mat<char>&,
        std::vector<int>&) const;
    %ignore getRawElementData;

    %pythoncode %{
    def getRawElementData(self, returnDomainIndices=False):
       if returnDomainIndices:
           return self._getRawElementDataWithElementIndices()
       else:
           return self._getRawElementDataNoElementIndices()
    %}

    %pythonappend indexSet %{
        val._parentGridView = self
    %}

    %pythonappend vtkWriter %{
        val._parentGridView = self
    %}

    %feature("compactdefaultargs") vtkWriter;
}

} // namespace Bempp

%include "grid/grid_view.hpp"

%clear arma::Mat<double>& vertices;
%clear arma::Mat<int>& elementCorners;
%clear arma::Mat<char>& auxData;
%clear std::vector<int>& domainIndices;
