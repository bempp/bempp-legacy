#ifndef BEMPP_PYTHON_SPACE_TYPES_H
#define BEMPP_PYTHON_SPACE_TYPES_H

#include<complex>
#include "bempp/common/shared_ptr.hpp"
#include "bempp/lib/grid.hpp"
#include "bempp/space/space.hpp"

<%
from space import dtypes
%>
namespace {
% for ctype in dtypes.itervalues():
%     if 'complex'  in ctype:
    typedef std::complex<${ctype[8:]}> ${ctype};
%     endif
% endfor
}
#endif
