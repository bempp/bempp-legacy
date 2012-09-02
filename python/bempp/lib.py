# Copyright (C) 2011-2012 by the BEM++ Authors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


import bempp.core as core

def checkType(dtype):
    dtypes={'float':'float64',
            'float32':'float32',
            'float64':'float64',
            'complex':'complex128',
            'complex64':'complex64',
            'complex128':'complex128'}
    if dtype in dtypes:
        return dtypes[dtype]
    else:
        raise ValueError('Data type does not exist.')

def promoteTypeToComplex(dtype):
    dtypes={'float':'complex128',
            'float32':'complex64',
            'float64':'complex128',
            'complex':'complex128',
            'complex64':'complex64',
            'complex128':'complex128'}
    if dtype in dtypes:
        return dtypes[dtype]
    else:
        raise ValueError('Data type does not exist.')

def _constructObjectTemplatedOnBasis(className, basisFunctionType, *args, **kwargs):
    fullName = className + "_" + checkType(basisFunctionType)
    try:
        class_ = getattr(core, fullName)
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

def _constructObjectTemplatedOnResult(className, resultType, *args, **kwargs):
    fullName = className + "_" + checkType(resultType)
    try:
        class_ = getattr(core, fullName)
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

def _constructObjectTemplatedOnValue(className, valueType,
                                    *args, **kwargs):
    fullName = className + "_" + checkType(valueType)
    try:
        class_ = getattr(core, fullName)
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

def _constructObjectTemplatedOnBasisAndResult(className,
                                              basisFunctionType, resultType,
                                              *args, **kwargs):
    # if basisFunctionType is None:
    #     if resultType is None:
    #         basisFunctionType = "float64"
    #         resultType = "float64"
    #     else:
    #         if resultType in ("float64", "complex128"):
    #             basisFunctionType = "float64"
    #         else:
    #             basisFunctionType = "float32"
    # else:
    #     if resultType is None:
    #         resultType = basisFunctionType

    fullName = (className + "_" +
                checkType(basisFunctionType) + "_" +
                checkType(resultType))
    try:
        class_ = getattr(core, fullName)
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

def _constructObjectTemplatedOnBasisKernelAndResult(
        className, basisFunctionType, kernelType, resultType, *args, **kwargs):
    fullName = (className + "_" +
                checkType(basisFunctionType) + "_" +
                checkType(kernelType) + "_" +
                checkType(resultType))
    try:
        class_ = getattr(core, fullName)
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

def createGridFactory():
    """Return a GridFactory object"""
    return core.GridFactory

def createNumericalQuadratureStrategy(basisFunctionType, resultType, accuracyOptions=None):
    """Construct a NumericalQuadratureStrategy object.

    *Arguments:*
        basisFunctionType (string)
            Type used to represent values of basis functions.

        resultType (string)
            Type used to represent values of boundary-element integrals.

        accuracyOptions (AccuracyOptions)
            Determines quadrature order. If set to None, default quadrature orders
            are used.

        The following combinations of basisFunctionType and resultType are allowed:

            basisFunctionType     resultType
            --------------------------------
            "float32"             "float32" or "complex64"
            "float64"             "float64" or "complex128"
            "complex64"           "complex64"
            "complex128"          "complex128"
    """
    basisFunctionType = checkType(basisFunctionType)
    resultType = checkType(resultType)
    if accuracyOptions is None:
        accuracyOptions = core.AccuracyOptions()
    name = 'numericalQuadratureStrategy'
    return _constructObjectTemplatedOnBasisAndResult(
        name, basisFunctionType, resultType, accuracyOptions)

def createAssemblyOptions(*args,**kwargs):
    """Return an AssemblyOptions object"""
    return core.AssemblyOptions(*args,**kwargs)

def createContext(factory, assemblyOptions):
    """Return an operator assembly context"""
    name = 'Context'
    return _constructObjectTemplatedOnBasisAndResult(
        name, factory.basisFunctionType(), factory.resultType(),
        factory, assemblyOptions)

def createPiecewiseConstantScalarSpace(context, grid):
    """Return a space of piecewise constant scalar functions"""
    name = 'piecewiseConstantScalarSpace'
    return _constructObjectTemplatedOnBasis(name, context.basisFunctionType(), grid)

def _constructOperator(className, context, domain, range, dualToRange):
    # determine basis function type
    basisFunctionType = domain.basisFunctionType()
    if (basisFunctionType != range.basisFunctionType() or
            basisFunctionType != dualToRange.basisFunctionType()):
        raise TypeError("BasisFunctionType of all spaces must be the same")

    # determine result type
    resultType = context.resultType()

    result = _constructObjectTemplatedOnBasisAndResult(
        className, basisFunctionType, resultType,
        context, domain, range, dualToRange)
    result._context = context
    result._domain = domain
    result._range = range
    result._dualToRange = dualToRange
    return result

def createPiecewiseLinearContinuousScalarSpace(context, grid):
    """Return space of piecewise linear continuous scalar functions"""
    name = 'piecewiseLinearContinuousScalarSpace'
    return _constructObjectTemplatedOnBasis(name, context.basisFunctionType(), grid)

def createLaplace3dSingleLayerBoundaryOperator(context, domain, range, dualToRange):
    """Return a single-layer-boundary operator for the Laplace equation in 3D."""
    return _constructOperator(
    "laplace3dSingleLayerBoundaryOperator", context, domain, range, dualToRange)

def createLaplace3dDoubleLayerBoundaryOperator(context, domain, range, dualToRange):
    """Return a double-layer-boundary operator for the Laplace equation in 3D."""
    return _constructOperator(
    "laplace3dDoubleLayerBoundaryOperator", context, domain, range, dualToRange)

def createLaplace3dAdjointDoubleLayerBoundaryOperator(context, domain, range, dualToRange):
    """Return an adjoint double-layer-boundary operator for the Laplace equation in 3D."""
    return _constructOperator(
    "laplace3dAdjointDoubleLayerBoundaryOperator", context, domain, range, dualToRange)

def createLaplace3dHypersingularBoundaryOperator(context, domain, range, dualToRange):
    """Return a hypersingular boundary operator for the Laplace equation in 3D."""
    return _constructOperator(
    "laplace3dHypersingularBoundaryOperator", context, domain, range, dualToRange)

def _constructLaplacePotentialOperator(className, context):
    basisFunctionType = context.basisFunctionType()
    resultType = context.resultType()
    result = _constructObjectTemplatedOnBasisAndResult(
        className, basisFunctionType, resultType)
    result._context = context
    return result

def createLaplace3dSingleLayerPotentialOperator(context):
    """Return a single-layer-potential operator for the Laplace equation in 3D."""
    return _constructLaplacePotentialOperator(
        "Laplace3dSingleLayerPotentialOperator", context)

def createLaplace3dDoubleLayerPotentialOperator(context):
    """Return a double-layer-potential operator for the Laplace equation in 3D."""
    return _constructLaplacePotentialOperator(
        "Laplace3dDoubleLayerPotentialOperator", context)

def _constructHelmholtzOperator(className, context, domain, range, dualToRange, waveNumber):
    basisFunctionType = context.basisFunctionType()
    if (basisFunctionType != domain.basisFunctionType() or
            basisFunctionType != range.basisFunctionType() or
            basisFunctionType != dualToRange.basisFunctionType()):
        raise TypeError("BasisFunctionType of context and all spaces must be the same")
    resultType = context.resultType()
    result = _constructObjectTemplatedOnBasis(
        className, basisFunctionType, context, domain, range, dualToRange, waveNumber)
    result._context = context
    result._domain = domain
    result._range = range
    result._dualToRange = dualToRange
    return result

def createHelmholtz3dSingleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Return a single-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "helmholtz3dSingleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def createHelmholtz3dDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Return a double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "helmholtz3dDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def createHelmholtz3dAdjointDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Return an adjoint double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "helmholtz3dAdjointDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def createHelmholtz3dHypersingularBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Return a hypersingular operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "helmholtz3dHypersingularBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def _constructHelmholtzPotentialOperator(className, context, waveNumber):
    basisFunctionType = context.basisFunctionType()
    resultType = context.resultType()
    result = _constructObjectTemplatedOnBasis(
        className, basisFunctionType, waveNumber)
    result._context = context
    return result

def createHelmholtz3dSingleLayerPotentialOperator(context, waveNumber):
    """Return a single-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzPotentialOperator(
        "helmholtz3dSingleLayerPotentialOperator", context, waveNumber)

def createHelmholtz3dDoubleLayerPotentialOperator(context, waveNumber):
    """Return a double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzPotentialOperator(
        "helmholtz3dDoubleLayerPotentialOperator", context, waveNumber)

def _constructModifiedHelmholtzOperator(className, context,
                                        domain, range, dualToRange, waveNumber):
    basisFunctionType = context.basisFunctionType()
    if (basisFunctionType != domain.basisFunctionType() or
            basisFunctionType != range.basisFunctionType() or
            basisFunctionType != dualToRange.basisFunctionType()):
        raise TypeError("BasisFunctionType of context and all spaces must be the same")
    resultType = context.resultType()

    waveNumberIsComplex = complex(waveNumber).imag != 0
    if waveNumberIsComplex and resultType in ("float32", "float64"):
        raise TypeError("Real result type given for a complex wave number")

    # determine kernelType
    if waveNumberIsComplex:
        kernelType = resultType
    else:
        if resultType in ("float32", "complex64"):
            kernelType = "float32"
        else:
            kernelType = "float64"

    # construct object
    result = _constructObjectTemplatedOnBasisKernelAndResult(
        className, basisFunctionType, kernelType, resultType,
        context, domain, range, dualToRange, waveNumber)
    result._context = context
    result._domain = domain
    result._range = range
    result._dualToRange = dualToRange
    return result

def createModifiedHelmholtz3dSingleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Return a single-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "modifiedHelmholtz3dSingleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def createModifiedHelmholtz3dDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Return a double-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "modifiedHelmholtz3dDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def createModifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Return an adjoint double-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "modifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def createModifiedHelmholtz3dHypersingularBoundaryOperator(
         context, domain, range, dualToRange, waveNumber):
     """Return a hypersingular operator for the modified Helmholtz equation in 3D."""
     return _constructModifiedHelmholtzOperator(
         "modifiedHelmholtz3dHypersingularBoundaryOperator", context, domain, range, dualToRange,
         waveNumber)


def createIdentityOperator(context, domain, range, dualToRange):
    """Return an identity operator."""
    return _constructOperator(
        "identityOperator", context, domain, range, dualToRange)

def __gridFunctionFromFunctor(
        functorType,
        context, space, dualSpace, function,
        argumentDimension, resultDimension):
    basisFunctionType = checkType(context.basisFunctionType())
    resultType = checkType(context.resultType())
    if (basisFunctionType != space.basisFunctionType() or
            basisFunctionType != dualSpace.basisFunctionType()):
        raise TypeError("BasisFunctionType of context, space and dualSpace must be the same")

    functor = _constructObjectTemplatedOnValue(
        "Python" + functorType,
        resultType, function, argumentDimension, resultDimension)
    result = _constructObjectTemplatedOnBasisAndResult(
        "gridFunctionFromPython" + functorType,
        basisFunctionType, resultType,
        context, space, dualSpace, functor)
    result._context = context
    result._space = space
    result._dualSpace = dualSpace
    return result

def gridFunctionFromSurfaceNormalDependentFunction(
        context, space, dualSpace, function,
        argumentDimension=3, resultDimension=1):
    return __gridFunctionFromFunctor(
        "SurfaceNormalDependentFunctor",
        context, space, dualSpace, function,
        argumentDimension, resultDimension)

def gridFunctionFromSurfaceNormalIndependentFunction(
        context, space, dualSpace, function,
        argumentDimension=3, resultDimension=1):
    return __gridFunctionFromFunctor(
        "SurfaceNormalIndependentFunctor",
        context, space, dualSpace, function,
        argumentDimension, resultDimension)

def createDefaultIterativeSolver(boundaryOperator,
                                 test_convergence="test_convergence_in_dual_to_range"):
    """Return the default iterative linear solver.

    This solver lets you solve the equation A f = g for the function
    f, with A being the boundary operator passed via the boundaryOperator
    argument and g a grid function supplied to the solve() method.
    """
    basisFunctionType = boundaryOperator.basisFunctionType()
    resultType = boundaryOperator.resultType()
    result = _constructObjectTemplatedOnBasisAndResult(
        "DefaultIterativeSolver", basisFunctionType, resultType,
        boundaryOperator,test_convergence)
    result._boundaryOperator = boundaryOperator
    return result

from core import defaultGmresParameterList
from core import defaultCgParameterList

def createAccuracyOptions():
    "Return an AccuracyOptions object"
    return core.AccuracyOptions()

def createAcaOptions():
    "Return an AcaOptions object"
    return core.AcaOptions()

def createEvaluationOptions():
    "Return and EvaluationOptions object"
    return core.EvaluationOptions()

def createBlockOperatorStructure(context):
    """Return an BlockedOperatorStructure object"""
    name = 'BlockedOperatorStructure'
    return _constructObjectTemplatedOnBasisAndResult(
        name, context.basisFunctionType(), context.resultType())

def createBlockedBoundaryOperator(context,structure):
    """Return a BlockedBoundaryOperator object"""
    name = 'BlockedBoundaryOperator'
    return _constructObjectTemplatedOnBasisAndResult(
        name, context.basisFunctionType(), context.resultType(),
        structure)

def createAcaApproximateLuInverse(operator,delta):
    """Return an AcaApproximateLuInverse Object"""
    name = 'createAcaApproximateLuInverse'
    return _constructObjectTemplatedOnValue(name,operator.valueType(),operator,delta)

def acaDiscreteOperatorToPreconditioner(operator,delta=1E-2):
    """Return an ACA Preconditioner"""
    name = 'acaDiscreteOperatorToPreconditioner'
    return _constructObjectTemplatedOnValue(name,operator.valueType(),operator,delta)

def acaBlockDiagonalPreconditioner(operators,deltas):
    """Return a block diagonal ACA Preconditioner"""
    name = 'acaBlockDiagonalPreconditioner'
    if len(operators)==0:
        raise TypeError("acaBlockDiagonalPreconditioner(): "
                        "Array 'operators' must not be empty.")
    typeName = operators[0].valueType()
    return _constructObjectTemplatedOnValue(name,typeName,operators,deltas)

def acaOperatorApproximateLuInverse(operator,delta):
    """Return the LU Approximate inverse of a DiscreteAcaOperator if it represents an H-Matrix"""
    name = 'acaOperatorApproximateLuInverse'
    return _constructObjectTemplatedOnValue(name,operator.valueType(),operator,delta)

def scaledAcaOperator(operator,multiplier):
    """Scale an H-Matrix by the factor given in multiplier"""
    name = 'scaledAcaOperator'
    return _constructObjectTemplatedOnValue(name,operator.valueType(),operator,multiplier)

def acaOperatorSum(op1,op2,eps,maximumRank):
    """Add two H-Matrices"""
    name = 'acaOperatorSum'
    if (op1.valueType() != op2.valueType()):
        raise TypeError("acaOperatorSum: ValueTypes of 'op1' and 'op2' do not match.")
    return _constructObjectTemplatedOnValue(name,op1.valueType(),op1,op2,eps,maximumRank)


