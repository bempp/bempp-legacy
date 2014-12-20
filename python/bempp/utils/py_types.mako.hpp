<%
from data_types import dtypes
%>\
#ifndef BEMPP_PYTHON_TYPES_H
#define BEMPP_PYTHON_TYPES_H

#include<complex>

namespace {
% for ctype in dtypes.values():
%     if 'complex'  in ctype:
    //! Typedef to avoid deep template nesting cython can't handle
    typedef std::complex<${ctype[8:]}> ${ctype};
%     endif
% endfor
}
#endif
