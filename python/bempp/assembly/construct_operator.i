%pythoncode %{

def _constructOperator(className, domain, range, dualToRange, resultType):
    # determine basis function type
    basisFunctionType = domain.basisFunctionType()
    if (basisFunctionType != range.basisFunctionType() or
            basisFunctionType != dualToRange.basisFunctionType()):
        raise TypeError("BasisFunctionType of all spaces must be the same")

    # determine result type
    if resultType is None:
        resultType = basisFunctionType
    resultType = checkType(resultType)

    result = constructObjectTemplatedOnBasisAndResult(
        className, basisFunctionType, resultType,
        domain, range, dualToRange)
    result._domain = domain
    result._range = range
    result._dualToRange = dualToRange
    return result

%}
