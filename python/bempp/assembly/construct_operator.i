%pythoncode %{

def _constructOperator(className, testSpace, trialSpace, resultType):
    basisFunctionType = testSpace.basisFunctionType()
    if resultType is None: resultType=basisFunctionType
    resultType=checkType(resultType)
    if (basisFunctionType != trialSpace.basisFunctionType()):
        raise TypeError("BasisFunctionType of testSpace must match that of trialSpace")
    return constructObjectTemplatedOnBasisAndResult(
        className, basisFunctionType, resultType,
        testSpace, trialSpace)

%}
