%{
#include "assembly/laplace_3d_single_layer_potential.hpp"
%}

// TODO
// %include "laplace_3d_single_layer_potential_docstrings.i"

%include "assembly/laplace_3d_single_layer_potential.hpp"

namespace Bempp
{
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dSingleLayerPotential);
}

%pythoncode %{

class Laplace3dSingleLayerPotential(Template2, ElementarySingularIntegralOperator):
    pass
    # TODO: docs

%}
