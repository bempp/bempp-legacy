%{
#include "grid/grid_segment.hpp"
%}

%include "std_set.i"
%include "std_vector.i"

#define shared_ptr boost::shared_ptr
%include "grid/grid_segment.hpp"
#undef shared_ptr
