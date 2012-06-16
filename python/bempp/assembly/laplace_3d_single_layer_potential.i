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

  def Laplace3dSingleLayerPotential(testSpace,trialSpace,BasisFunctionType=None,ResultType=None):
      if BasisFunctionType is None: BasisFunctionType=testSpace._basisFunctionType
      if ResultType is None: ResultType=BasisFunctionType
      dtype1=checkType(BasisFunctionType)
      dtype2=checkType(ResultType)
      name='Laplace3dSingleLayerPotential'+'_'+dtype1+'_'+dtype2
      fullName=name+'_'+dtype1+'_'+dtype2
      keys=globals()
      if not name in keys: 
          raise Exception('Class '+name+' with BasisFunctionType '+dtype1+' and ResultType '+dtype2+' does not exist.')
      return keys[fullName](testSpace,trialSpace)


%}
