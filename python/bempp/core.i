// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// Construct a Swig module
%module core
%{
#define SWIG_FILE_WITH_INIT

#include <dune/common/exceptions.hh>
#include <complex>
%}

// Import docstring macros
%include "docstrings.i"

// Define commonly used typemaps
%include "armadillo.i"
%include "auto_ptr.i"
%include "exception.i"
%include "std_string.i"

// Useful Python tools
%include "py_defs.i"

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
%inline %{
    namespace Bempp
    {
        typedef double ctype;
    }
%}
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

 // Space
 %include "space/space.i"
 %include "space/scalar_space.i"
 %include "space/piecewise_constant_scalar_space.i"
