%pythoncode %{

def _constructOperator(className, testSpace, trialSpace, resultType):
    basisFunctionType = testSpace.basisFunctionType()
    if resultType is None: resultType=basisFunctionType
    resultType=checkType(resultType)
    if (basisFunctionType != trialSpace.basisFunctionType()):
        raise TypeError("BasisFunctionType of testSpace must match that of trialSpace")
    result = constructObjectTemplatedOnBasisAndResult(
        className, basisFunctionType, resultType,
        testSpace, trialSpace)
    result._testSpace = testSpace
    result._trialSpace = trialSpace
    return result

%}
