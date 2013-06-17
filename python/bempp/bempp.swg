%module(package="bempp") core
// %include "core.i"
%{
#include <dune/common/exceptions.hh>
%}
%include "exception.i"
%import "core.i"
