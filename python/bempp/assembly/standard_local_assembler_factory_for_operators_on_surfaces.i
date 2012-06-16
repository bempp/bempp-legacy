%{
#include "assembly/standard_local_assembler_factory_for_operators_on_surfaces.hpp"
  %}


%include "assembly/standard_local_assembler_factory_for_operators_on_surfaces.hpp"


namespace Bempp
{
BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(StandardLocalAssemblerFactoryForOperatorsOnSurfaces);


}





%pythoncode %{

  def standardLocalAssemblerFactoryForOperatorsOnSurfaces(accuracyOptions=None,
      basisFunctionType='float64',resultType='float64'):
      if basisFunctionType is not None: dtype1=checkType(basisFunctionType)
      if resultType is not None: dtype2=checkType(resultType)
      if accuracyOptions is None: accuracyOptions=AccuracyOptions()
      name='StandardLocalAssemblerFactoryForOperatorsOnSurfaces'
      return constructObjectTemplatedOnBasisAndResult(name, basisFunctionType, resultType,accuracyOptions)
		      
	  %}
