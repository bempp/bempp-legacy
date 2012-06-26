%{
#include "assembly/default_local_assembler_factory_for_operators_on_surfaces.hpp"
%}

%include "assembly/default_local_assembler_factory_for_operators_on_surfaces.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
DefaultLocalAssemblerFactoryForOperatorsOnSurfaces);
}

%pythoncode %{

def defaultLocalAssemblerFactoryForOperatorsOnSurfaces(accuracyOptions=None,
        basisFunctionType='float64',resultType='float64'):
    if basisFunctionType is not None:
        basisFunctionType = checkType(basisFunctionType)
    if resultType is not None:
        resultType = checkType(resultType)
    if accuracyOptions is None:
        accuracyOptions = AccuracyOptions()
    name = 'DefaultLocalAssemblerFactoryForOperatorsOnSurfaces'
    return constructObjectTemplatedOnBasisAndResult(
        name, basisFunctionType, resultType, accuracyOptions)

%}
