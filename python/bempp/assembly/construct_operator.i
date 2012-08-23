%pythoncode %{

def _constructOperator(className, context, domain, range, dualToRange):
    # determine basis function type
    basisFunctionType = domain.basisFunctionType()
    if (basisFunctionType != range.basisFunctionType() or
            basisFunctionType != dualToRange.basisFunctionType()):
        raise TypeError("BasisFunctionType of all spaces must be the same")

    # determine result type
    resultType = context.resultType()

    result = constructObjectTemplatedOnBasisAndResult(
        className, basisFunctionType, resultType,
        context, domain, range, dualToRange)
    result._context = context
    result._domain = domain
    result._range = range
    result._dualToRange = dualToRange
    return result

%}
