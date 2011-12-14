// Construct a Swig module
%module bempp
%{
#define SWIG_FILE_WITH_INIT

#include <dune/common/exceptions.hh>
%}

// Import docstring macros
%include "docstrings.i"

// Define commonly used typemaps
%include "armadillo.i"
%include "auto_ptr.i"
%include "exception.i"
%include "std_string.i"

// Setup a handler for C++ exceptions
%exception {
    try {
        $action 
    }
    catch (const Dune::Exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what().c_str());
    }
    catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}

// Declare all necessary auto_ptr typemaps
AUTO_PTR_TYPEMAPS(Bempp::Grid)
AUTO_PTR_TYPEMAPS(Bempp::GridView)
AUTO_PTR_TYPEMAPS(Bempp::EntityPointer<0>)
AUTO_PTR_TYPEMAPS(Bempp::EntityPointer<1>)
AUTO_PTR_TYPEMAPS(Bempp::EntityPointer<2>)
AUTO_PTR_TYPEMAPS(Bempp::EntityPointer<3>)
AUTO_PTR_TYPEMAPS(Bempp::EntityIterator<0>)
AUTO_PTR_TYPEMAPS(Bempp::EntityIterator<1>)
AUTO_PTR_TYPEMAPS(Bempp::EntityIterator<2>)
AUTO_PTR_TYPEMAPS(Bempp::EntityIterator<3>)
AUTO_PTR_TYPEMAPS(Bempp::VtkWriter)

// Make commonly used typedefs known to Swig
%include "grid/common.hpp"
%include "grid/geometry_type.hpp"

// Wrap Bempp components

// Grid
%include "grid/geometry.i"  
%include "grid/geometry_type.i"  
%include "grid/entity.i"  
%include "grid/entity_pointer.i"  
%include "grid/entity_iterator.i"  
%include "grid/grid.i"  
%include "grid/id_set.i"
%include "grid/grid_view.i"  
%include "grid/index_set.i"
%include "grid/vtk_writer.i" 
%include "grid/grid_factory.i" 
