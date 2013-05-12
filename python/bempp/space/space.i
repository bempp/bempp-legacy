%{
#include "space/space.hpp"
#include <complex>
%}

%shared_ptr(Bempp::Space<float>);
%shared_ptr(Bempp::Space<double>);
%shared_ptr(Bempp::Space<std::complex<float> >);
%shared_ptr(Bempp::Space<std::complex<double> >);

// TODO
// %include "space_docstrings.i"

namespace Bempp
{

template<typename BasisFunctionType> class Space;

%extend Space
{

// this function is only for internal use
%ignore basis;

// this function is only for internal use
%ignore shapeFunctionValue;

// to be wrapped later...
%ignore setElementVariant;

// to be wrapped later...
%ignore elementVariant;

%ignore getGlobalDofs;

%ignore global2localDofs;
%ignore flatLocal2localDofs;

// these functions are only for internal use
%ignore getGlobalDofPositions;
%ignore getFlatLocalDofPositions;
%ignore getGlobalDofBoundingBoxes;
%ignore getFlatLocalDofBoundingBoxes;
%ignore getGlobalDofNormals;
%ignore getFlatLocalDofNormals;

%ignore dumpClusterIds;
%ignore dumpClusterIdsEx;
}

%define BEMPP_EXTEND_SPACE(BASIS, PYBASIS)
    %extend Space< BASIS >
    {
        %pythonprepend assignDofs %{
            print ("HINT: It is not necessary to call Space.assignDofs() any more.\n"
                   "DOFs are now automatically assigned as soon as a space object is "
                   "constructed.")
        %}
    }
%enddef
BEMPP_ITERATE_OVER_BASIS_TYPES(BEMPP_EXTEND_SPACE);

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS(Space);

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "space/space.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(Space);
}

