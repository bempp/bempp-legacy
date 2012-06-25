%{
#include "assembly/standard_local_assembler_factory_for_operators_on_surfaces.hpp"
%}

%include "assembly/standard_local_assembler_factory_for_operators_on_surfaces.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
StandardLocalAssemblerFactoryForOperatorsOnSurfaces);
}

%pythoncode %{

def standardLocalAssemblerFactoryForOperatorsOnSurfaces(accuracyOptions=None,
        basisFunctionType='float64',resultType='float64'):
    if basisFunctionType is not None:
        basisFunctionType = checkType(basisFunctionType)
    if resultType is not None:
        resultType = checkType(resultType)
    if accuracyOptions is None:
        accuracyOptions = AccuracyOptions()
    name = 'StandardLocalAssemblerFactoryForOperatorsOnSurfaces'
    return constructObjectTemplatedOnBasisAndResult(
        name, basisFunctionType, resultType, accuracyOptions)

%}
