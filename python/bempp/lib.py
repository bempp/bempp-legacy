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
import numpy

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

def _constructObjectTemplatedOnBasis(module, className, basisFunctionType,
                                     *args, **kwargs):
    fullName = className + "_" + checkType(basisFunctionType)
    try:
        class_ = getattr(module, fullName)
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

def _constructObjectTemplatedOnResult(module, className, resultType,
                                      *args, **kwargs):
    fullName = className + "_" + checkType(resultType)
    try:
        class_ = getattr(module, fullName)
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

def _constructObjectTemplatedOnValue(module, className, valueType,
                                     *args, **kwargs):
    fullName = className + "_" + checkType(valueType)
    try:
        class_ = getattr(module, fullName)
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

def _constructObjectTemplatedOnBasisAndResult(module, className,
                                              basisFunctionType, resultType,
                                              *args, **kwargs):
    fullName = (className + "_" +
                checkType(basisFunctionType) + "_" +
                checkType(resultType))
    try:
        class_ = getattr(module, fullName)
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

def _constructObjectTemplatedOnBasisKernelAndResult(
        module, className, basisFunctionType, kernelType, resultType,
        *args, **kwargs):
    fullName = (className + "_" +
                checkType(basisFunctionType) + "_" +
                checkType(kernelType) + "_" +
                checkType(resultType))
    try:
        class_ = getattr(module, fullName)
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

def createGridFactory():
    """Return a GridFactory object"""
    return core.GridFactory

def createNumericalQuadratureStrategy(basisFunctionType, resultType, accuracyOptions=None):
    """
    Create and return a NumericalQuadratureStrategy object.

    A quadrature strategy provides functions constructing local assemblers used
    to discretize boundary operators and user-defined functions. A particular
    quadrature strategy determines how the integrals involved in this
    discretization are evaluated.

    The local assemblers constructed by this class use numerical quadrature to
    evaluate the necessary integrals. Singular integrals are transformed into
    regular ones as described in S. Sauter, Ch. Schwab, "Boundary Element
    Methods" (2010). Quadrature accuracy can be influenced by the
    'accuracyOptions' parameter.

    *Parameters:*
       - basisFunctionType (string)
            Type used to represent the values of the (components of the) basis
            functions into which arguments of operators discretized with this
            strategy will be expanded.

       - resultType (string)
            Type used to represent the values of integrals.

       - accuracyOptions (AccuracyOptions or AccuracyOptionsEx)
            Determines quadrature orders used to approximate different types of
            integrals. If set to None, default quadrature orders are used.

    The following combinations of basisFunctionType and resultType are allowed:

    =================     =========================
    basisFunctionType     resultType
    =================     =========================
    "float32"             "float32" or "complex64"
    "float64"             "float64" or "complex128"
    "complex64"           "complex64"
    "complex128"          "complex128"
    =================     =========================

    Typically, you should choose basisFunctionType = "float64" and resultType =
    "float64" for calculations with real-valued operators (such as those related
    to the Laplace equation) or basisFunctionType = "float64" and resultType =
    "complex128" for calculations with complex-valued operators (such as those
    related to the Helmlholtz equation).
    """
    basisFunctionType = checkType(basisFunctionType)
    resultType = checkType(resultType)
    if accuracyOptions is None:
        accuracyOptions = core.AccuracyOptions()
    name = 'numericalQuadratureStrategy'
    return _constructObjectTemplatedOnBasisAndResult(
        core, name, basisFunctionType, resultType, accuracyOptions)

def createAssemblyOptions():
    """Create and return an AssemblyOptions object with default settings."""
    return core.AssemblyOptions()

def createContext(factory, assemblyOptions):
    """
    Create and return a Context object.

    A Context determines the mechanics of the assembly of weak forms and the
    evaluation of potentials.

    *Parameters:*
       - quadStrategy (QuadratureStrategy)
            Quadrature strategy to be used for the calculation of integrals
            occurring e.g. in the weak forms of boundary operators or in the
            definition of potential operators.

       - assemblyOptions (AssemblyOptions)
            Further options influencing the weak-form assembly process.

    *Returns* a newly constructed Context_BasisFunctionType_ResultType object,
    with BasisFunctionType and ResultType determined automatically from the
    quadStrategy argument and equal to either float32, float64, complex64 or
    complex128.
    """
    name = 'Context'
    return _constructObjectTemplatedOnBasisAndResult(
        core, name, factory.basisFunctionType(), factory.resultType(),
        factory, assemblyOptions)

def createPiecewiseConstantScalarSpace(context, grid, segment=None):
    """
    Create and return a space of scalar functions defined on a grid (or its
    segment) and constant on each element of this grid.

    *Parameters:*
       - context (Context)
            A Context object that will determine the type used to represent the
            values of the basis functions of the newly constructed space.
       - grid (Grid)
            Grid on which the functions from the newly constructed space will be
            defined.
       - segment (GridSegment)
            (Optional) Segment of the grid on which the space should be defined.
            If set to None (default), the whole grid will be used.

    *Returns* a newly constructed Space_BasisFunctionType object, with
    BasisFunctionType determined automatically from the context argument and
    equal to either float32, float64, complex64 or complex128.
    """
    name = 'piecewiseConstantScalarSpace'
    return _constructObjectTemplatedOnBasis(
        core, name, context.basisFunctionType(), grid, segment)

def createPiecewiseConstantDualMeshScalarSpaceBarycentric(context, grid, segment=None):
    """
    Create and return a space of scalar functions defined on the dual of a grid (or its
    segment) and constant on each element of the dual grid.

    *Parameters:*
       - context (Context)
            A Context object that will determine the type used to represent the
            values of the basis functions of the newly constructed space.
       - grid (Grid)
            Grid on which the functions from the newly constructed space will be
            defined. The grid must have a barycentric refinement level.
       - segment (GridSegment)
            (Optional) Segment of the grid on which the space should be defined.
            If set to None (default), the whole grid will be used.

       *Returns* a newly constructed Space_BasisFunctionType object, with
       BasisFunctionType determined automatically from the context argument and
       equal to either float32, float64, complex64 or complex128.
    """
    name = 'piecewiseConstantDualMeshScalarSpaceBarycentric'
    return _constructObjectTemplatedOnBasis(
          core, name, context.basisFunctionType(), grid, segment)


def createPiecewiseLinearContinuousScalarSpace(context, grid, segment=None,
                                               strictlyOnSegment=False):
    """
    Create and return a space of globally continuous scalar functions defined
    on a grid (or its segment) and linear on each element of this grid.

    *Parameters:*
       - context (Context)
            A Context object that will determine the type used to represent the
            values of the basis functions of the newly constructed space.
       - grid (Grid)
            Grid on which the functions from the newly constructed space will be
            defined.
       - segment (GridSegment)
            (Optional) Segment of the grid on which the space should be defined.
            If set to None (default), the whole grid will be used.

    *Returns* a newly constructed Space_BasisFunctionType object, with
    BasisFunctionType determined automatically from the context argument and
    equal to either float32, float64, complex64 or complex128.
    """
    name = 'piecewiseLinearContinuousScalarSpace'
    return _constructObjectTemplatedOnBasis(
        core, name, context.basisFunctionType(), grid, segment,
        strictlyOnSegment)

def createPiecewiseLinearContinuousScalarSpaceBarycentric(context, grid, segment=None,
                                                          strictlyOnSegment=False):
    """
    Create and return a space of globally continuous scalar functions defined
    on grid (or its segment) and linear on each element of this grid. It produces the
    same global dofs as *createPiecewiseLinearContinuousScalarSpace*, but the local
    dofs live on a barycentric refinement of the original grid. Hence, this space is
    compatible with other spaces defined over barycentric refinements.

    *Parameters:*
       - context (Context)
            A Context object that will determine the type used to represent the
            values of the basis functions of the newly constructed space.
       - grid (Grid)
            Grid on which the functions from the newly constructed space will be
            defined.
       - segment (GridSegment)
            (Optional) Segment of the grid on which the space should be defined.
            If set to None (default), the whole grid will be used.

    *Returns* a newly constructed Space_BasisFunctionType object, with
    BasisFunctionType determined automatically from the context argument and
    equal to either float32, float64, complex64 or complex128.
    """
    name = 'piecewiseLinearContinuousScalarSpaceBarycentric'
    return _constructObjectTemplatedOnBasis(
          core, name, context.basisFunctionType(), grid, segment,
          strictlyOnSegment)


def createPiecewiseLinearDiscontinuousScalarSpace(
    context, grid, segment=None, strictlyOnSegment=False):
    """
    Create and return a space of scalar functions defined on a grid (or its
    segment) and linear on each element of this grid (but not forced to be
    continuous at element boundaries).

    *Parameters:*
       - context (Context)
            A Context object that will determine the type used to represent the
            values of the basis functions of the newly constructed space.
       - grid (Grid)
            Grid on which the functions from the newly constructed space will be
            defined.
       - segment (GridSegment)
            (Optional) Segment of the grid on which the space should be defined.
            If set to None (default), the whole grid will be used.
       - strictlyOnSegment (bool)
            (Optional) If set to False (default), the space will include all
            basis functions associated with vertices belonging to the chosen
            segment, regardless of whether the elements on which these functions
            are defined belong themselves to the segment. As a result, the
            resulting space will be (in the mathematical sense) a superset of a
            PiecewiseLinearContinuousScalarSpace defined on the same segment. If
            set to True, the space will only include basis functions defined on
            elements belonging to the chosen segment.

    *Returns* a newly constructed Space_BasisFunctionType object, with
    BasisFunctionType determined automatically from the context argument and
    equal to either float32, float64, complex64 or complex128.
    """
    name = 'piecewiseLinearDiscontinuousScalarSpace'
    return _constructObjectTemplatedOnBasis(
        core, name, context.basisFunctionType(), grid,
        segment, strictlyOnSegment)

def createPiecewisePolynomialContinuousScalarSpace(
        context, grid, polynomialOrder, segment=None,
        strictlyOnSegment=False):
    """
    Create and return a space of globally continuous scalar functions defined on
    a grid (or its segment) and having a polynomial representation of a given
    order on each element of this grid.

    *Parameters:*
       - context (Context)
            A Context object that will determine the type used to represent the
            values of the basis functions of the newly constructed space.
       - grid (Grid)
            Grid on which the functions from the newly constructed space will be
            defined.
       - polynomialOrder (int)
            Order of the polynomial basis defined on each element.
       - segment (GridSegment)
            (Optional) Segment of the grid on which the space should be defined.
            If set to None (default), the whole grid will be used.

    *Returns* a newly constructed Space_BasisFunctionType object, with
    BasisFunctionType determined automatically from the context argument and
    equal to either float32, float64, complex64 or complex128.
    """
    name = 'piecewisePolynomialContinuousScalarSpace'
    return _constructObjectTemplatedOnBasis(
        core, name, context.basisFunctionType(), grid, polynomialOrder,
        segment, strictlyOnSegment)

def createPiecewisePolynomialDiscontinuousScalarSpace(
        context, grid, polynomialOrder, segment=None,
        requireReferencePointOnSegment=True,
        requireElementOnSegment=False):
    """
    Create and return a space of scalar functions defined on a grid (or its
    segment) and having a polynomial representation of a given order on each
    element of this grid (but not forced to be continuous at element
    boundaries).

    *Parameters:*
       - context (Context)
            A Context object that will determine the type used to represent the
            values of the basis functions of the newly constructed space.
       - grid (Grid)
            Grid on which the functions from the newly constructed space will be
            defined.
       - polynomialOrder (int)
            Order of the polynomial basis defined on each element.
       - segment (GridSegment)
            (Optional) Segment of the grid on which the space should be defined.
            If set to None (default), the whole grid will be used.
       - strictlyOnSegment (bool)
            (Optional) If set to False (default), the space will include those
            and only those basis functions that are necessary to make the newly
            constructed space a superset (in the mathematical sense) of a
            PiecewisePolynomialContinuousScalarSpace defined on the chosen
            segment. Otherwise the space will include only the functions defined
            on elements belonging to the segment.

    *Returns* a newly constructed Space_BasisFunctionType object, with
    BasisFunctionType determined automatically from the context argument and
    equal to either float32, float64, complex64 or complex128.
    """
    name = 'piecewisePolynomialDiscontinuousScalarSpace'
    dofMode = 0
    if requireReferencePointOnSegment:
        dofMode |= 4
    if requireElementOnSegment:
        dofMode |= 2
    return _constructObjectTemplatedOnBasis(
        core, name, context.basisFunctionType(), grid, polynomialOrder,
        segment, dofMode)

def createUnitScalarSpace(context, grid):
    """
    Create and return a space consisting of the single function equal to 1
    on all elements of the grid.

    *Parameters:*
       - context (Context)
            A Context object that will determine the type used to represent the
            values of the basis functions of the newly constructed space.
       - grid (Grid)
            Grid on which the function from the newly constructed space will be
            defined.

    *Returns* a newly constructed Space_BasisFunctionType object, with
    BasisFunctionType determined automatically from the context argument and
    equal to either float32, float64, complex64 or complex128.
    """
    return _constructObjectTemplatedOnBasis(
        core, 'unitScalarSpace', context.basisFunctionType(), grid)

def createRaviartThomas0VectorSpace(context, grid, segment=None,
                                    putDofsOnBoundaries=False,
                                    requireEdgeOnSegment=True,
                                    requireElementOnSegment=False):
    """
    Create and return a space of lowest order Raviart-Thomas vector functions
    with normal components continuous on boundaries between elements.

    *Parameters:*
       - context (Context)
            A Context object that will determine the type used to represent the
            values of the basis functions of the newly constructed space.
       - grid (Grid)
            Grid on which the functions from the newly constructed space will be
            defined.
       - segment (GridSegment)
            (Optional) Segment of the grid on which the space should be defined.
            If set to None (default), the whole grid will be used.
       - putDofsOnBoundaries (bool)
            (Optional) If set to False (default), degrees of freedom will not be
            placed on edges lying on boundaries of the grid. This is usually the
            desired behaviour for simulations of open perfectly conducting
            surfaces (sheets). If set to True, degrees of freedom will be placed
            on all edges belonging to the chosen segment of the grid.

    *Returns* a newly constructed Space_BasisFunctionType object, with
    BasisFunctionType determined automatically from the context argument and
    equal to either float32, float64, complex64 or complex128.
    """
    name = 'raviartThomas0VectorSpace'
    dofMode = 0
    if requireEdgeOnSegment:
        dofMode |= 1
    if requireElementOnSegment:
        dofMode |= 2
    return _constructObjectTemplatedOnBasis(
        core, name, context.basisFunctionType(), grid, segment,
        putDofsOnBoundaries, dofMode)

def _constructOperator(className, context, domain, range, dualToRange, label=None):
    # determine basis function type
    basisFunctionType = domain.basisFunctionType()
    if (basisFunctionType != range.basisFunctionType() or
            basisFunctionType != dualToRange.basisFunctionType()):
        raise TypeError("BasisFunctionType of all spaces must be the same")

    # determine result type
    resultType = context.resultType()

    if label:
        result = _constructObjectTemplatedOnBasisAndResult(
            core, className, basisFunctionType, resultType,
            context, domain, range, dualToRange, label)
    else:
        result = _constructObjectTemplatedOnBasisAndResult(
            core, className, basisFunctionType, resultType,
            context, domain, range, dualToRange)
    result._context = context
    result._domain = domain
    result._range = range
    result._dualToRange = dualToRange
    return result

# determine the type used to represent the values of the basis functions into
# which functions acted upon by the operator will be expanded, and the type used
# to represent the values of the functions produced by this operator.

def createLaplace3dSingleLayerBoundaryOperator(
        context, domain, range, dualToRange, label=None):
    """
    Create and return a single-layer-potential boundary operator for the
    Laplace equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructOperator(
        "laplace3dSingleLayerBoundaryOperator",
        context, domain, range, dualToRange, label)

def createLaplace3dDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, label=None):
    """
    Create and return a double-layer-potential boundary operator for the
    Laplace equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructOperator(
        "laplace3dDoubleLayerBoundaryOperator",
        context, domain, range, dualToRange, label)

def createLaplace3dAdjointDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, label=None):
    """
    Create and return an adjoint double-layer-potential boundary operator for
    the Laplace equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructOperator(
        "laplace3dAdjointDoubleLayerBoundaryOperator",
        context, domain, range, dualToRange, label)

def createLaplace3dHypersingularBoundaryOperator(
        context, domain, range, dualToRange, label=None):
    """
    Create and return a hypersingular boundary operator for the
    Laplace equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructOperator(
        "laplace3dHypersingularBoundaryOperator",
        context, domain, range, dualToRange, label)

def _constructLaplacePotentialOperator(className, context):
    basisFunctionType = context.basisFunctionType()
    resultType = context.resultType()
    result = _constructObjectTemplatedOnBasisAndResult(
        core, className, basisFunctionType, resultType)
    result._context = context
    return result

def createLaplace3dSingleLayerPotentialOperator(context):
    """
    Create and return a single-layer potential operator for the Laplace
    equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object used to control the evaluation of integrals
            occurring in the definition of the potential operator.

    *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.

    Note about BEM++ terminology: a *potential operator* acts on functions
    defined on a surface S and produces functions defined at any point of the
    space surrounding S, but not necessarily on S itself. In contrast, a
    *boundary operator* acts on on functions defined on a surface S and produces
    functions defined on the same surface S.
    """
    return _constructLaplacePotentialOperator(
        "laplace3dSingleLayerPotentialOperator", context)

def createLaplace3dDoubleLayerPotentialOperator(context):
    """
    Create and return a double-layer potential operator for the Laplace
    equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object used to control the evaluation of integrals
            occurring in the definition of the potential operator.

    *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.

    Note about BEM++ terminology: a *potential operator* acts on functions
    defined on a surface S and produces functions defined at any point of the
    space surrounding S, but not necessarily on S itself. In contrast, a
    *boundary operator* acts on on functions defined on a surface S and produces
    functions defined on the same surface S.
    """
    return _constructLaplacePotentialOperator(
        "laplace3dDoubleLayerPotentialOperator", context)

def _constructHelmholtzOperator(
        className, context, domain, range, dualToRange, waveNumber,
        label, useInterpolation, interpPtsPerWavelength):
    basisFunctionType = context.basisFunctionType()
    if (basisFunctionType != domain.basisFunctionType() or
            basisFunctionType != range.basisFunctionType() or
            basisFunctionType != dualToRange.basisFunctionType()):
        raise TypeError("BasisFunctionType of 'context' and all spaces "
                        "must be the same")
    resultType = context.resultType()
    if resultType != promoteTypeToComplex(basisFunctionType):
        raise TypeError("ResultType of 'context' must be the class of complex "
                        "numbers of the same precision as its "
                        "BasisFunctionType")
    if not label:
        label = ""
    symmetry = 0
    result = _constructObjectTemplatedOnBasis(
        core, className, basisFunctionType, context, domain, range, dualToRange,
        waveNumber, label, symmetry, useInterpolation, interpPtsPerWavelength)
    result._context = context
    result._domain = domain
    result._range = range
    result._dualToRange = dualToRange
    return result

def createHelmholtz3dSingleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber, label=None,
        useInterpolation=False, interpPtsPerWavelength=5000):
    """
    Create and return a single-layer-potential boundary operator for
    the Helmholtz equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the Helmholtz equation
                nabla^2 u + k^2 u = 0.
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.
       - useInterpolation (bool)
            If set to False (default), the standard exp() function will be used
            to evaluate the exponential factor occurring in the kernel. If set
            to True, the exponential factor will be evaluated by piecewise-cubic
            interpolation of values calculated in advance on a regular
            grid. This normally speeds up calculations, but might result in a
            loss of accuracy. This is an experimental feature: use it at your
            own risk.
       - interpPtsPerWavelength (int)
            If useInterpolation is set to True, this parameter determines the
            number of points per "effective wavelength" (defined as 2 pi /
            abs(waveNumber)) used to construct the interpolation grid. The
            default value (5000) is normally enough to reduce the relative or
            absolute error, *whichever is smaller*, below 100 * machine
            precision. If useInterpolation is set to False, this parameter is
            ignored.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructHelmholtzOperator(
        "helmholtz3dSingleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber,
        label, useInterpolation, interpPtsPerWavelength)

def createHelmholtz3dDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber, label=None,
        useInterpolation=False, interpPtsPerWavelength=5000):
    """
    Create and return a double-layer-potential boundary operator for
    the Helmholtz equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the Helmholtz equation
                nabla^2 u + k^2 u = 0.
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.
       - useInterpolation (bool)
            If set to False (default), the standard exp() function will be used
            to evaluate the exponential factor occurring in the kernel. If set
            to True, the exponential factor will be evaluated by piecewise-cubic
            interpolation of values calculated in advance on a regular
            grid. This normally speeds up calculations, but might result in a
            loss of accuracy. This is an experimental feature: use it at your
            own risk.
       - interpPtsPerWavelength (int)
            If useInterpolation is set to True, this parameter determines the
            number of points per "effective wavelength" (defined as 2 pi /
            abs(waveNumber)) used to construct the interpolation grid. The
            default value (5000) is normally enough to reduce the relative or
            absolute error, *whichever is smaller*, below 100 * machine
            precision. If useInterpolation is set to False, this parameter is
            ignored.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructHelmholtzOperator(
        "helmholtz3dDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber,
        label, useInterpolation, interpPtsPerWavelength)

def createHelmholtz3dAdjointDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber, label=None,
        useInterpolation=False, interpPtsPerWavelength=5000):
    """
    Create and return an adjoint double-layer-potential boundary operator for
    the Helmholtz equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the Helmholtz equation
                nabla^2 u + k^2 u = 0.
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.
       - useInterpolation (bool)
            If set to False (default), the standard exp() function will be used
            to evaluate the exponential factor occurring in the kernel. If set
            to True, the exponential factor will be evaluated by piecewise-cubic
            interpolation of values calculated in advance on a regular
            grid. This normally speeds up calculations, but might result in a
            loss of accuracy. This is an experimental feature: use it at your
            own risk.
       - interpPtsPerWavelength (int)
            If useInterpolation is set to True, this parameter determines the
            number of points per "effective wavelength" (defined as 2 pi /
            abs(waveNumber)) used to construct the interpolation grid. The
            default value (5000) is normally enough to reduce the relative or
            absolute error, *whichever is smaller*, below 100 * machine
            precision. If useInterpolation is set to False, this parameter is
            ignored.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructHelmholtzOperator(
        "helmholtz3dAdjointDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber,
        label, useInterpolation, interpPtsPerWavelength)

def createHelmholtz3dHypersingularBoundaryOperator(
        context, domain, range, dualToRange, waveNumber, label=None,
        useInterpolation=False, interpPtsPerWavelength=5000):
    """
    Create and return a hypersingular boundary operator for
    the Helmholtz equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the Helmholtz equation
                nabla^2 u + k^2 u = 0.
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.
       - useInterpolation (bool)
            If set to False (default), the standard exp() function will be used
            to evaluate the exponential factor occurring in the kernel. If set
            to True, the exponential factor will be evaluated by piecewise-cubic
            interpolation of values calculated in advance on a regular
            grid. This normally speeds up calculations, but might result in a
            loss of accuracy. This is an experimental feature: use it at your
            own risk.
       - interpPtsPerWavelength (int)
            If useInterpolation is set to True, this parameter determines the
            number of points per "effective wavelength" (defined as 2 pi /
            abs(waveNumber)) used to construct the interpolation grid. The
            default value (5000) is normally enough to reduce the relative or
            absolute error, *whichever is smaller*, below 100 * machine
            precision. If useInterpolation is set to False, this parameter is
            ignored.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructHelmholtzOperator(
        "helmholtz3dHypersingularBoundaryOperator", context,
        domain, range, dualToRange, waveNumber,
        label, useInterpolation, interpPtsPerWavelength)

def _constructHelmholtzPotentialOperator(className, context, waveNumber):
    basisFunctionType = context.basisFunctionType()
    resultType = context.resultType()
    result = _constructObjectTemplatedOnBasis(
        core, className, basisFunctionType, waveNumber)
    result._context = context
    return result

def createHelmholtz3dSingleLayerPotentialOperator(context, waveNumber):
    """
    Create and return a single-layer potential operator for the Helmholtz
    equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object used to control the evaluation of integrals
            occurring in the definition of the potential operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the Helmholtz equation
                nabla^2 u + k^2 u = 0.

    *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.

    Note about BEM++ terminology: a *potential operator* acts on functions
    defined on a surface S and produces functions defined at any point of the
    space surrounding S, but not necessarily on S itself. In contrast, a
    *boundary operator* acts on on functions defined on a surface S and produces
    functions defined on the same surface S.
    """
    return _constructHelmholtzPotentialOperator(
        "helmholtz3dSingleLayerPotentialOperator", context, waveNumber)

def createHelmholtz3dDoubleLayerPotentialOperator(context, waveNumber):
    """
    Create and return a double-layer potential operator for the Helmholtz
    equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object used to control the evaluation of integrals
            occurring in the definition of the potential operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the Helmholtz equation
                nabla^2 u + k^2 u = 0.

    *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.

    Note about BEM++ terminology: a *potential operator* acts on functions
    defined on a surface S and produces functions defined at any point of the
    space surrounding S, but not necessarily on S itself. In contrast, a
    *boundary operator* acts on on functions defined on a surface S and produces
    functions defined on the same surface S.
    """
    return _constructHelmholtzPotentialOperator(
        "helmholtz3dDoubleLayerPotentialOperator", context, waveNumber)

def createHelmholtz3dFarFieldSingleLayerPotentialOperator(context, waveNumber):
    """
    Create and return a potential operator used to calculate part of the far-field
    pattern for a radiating solution of the Helmholtz equation in 3D.

    See the C++ documentation for mathematical background.

    *Parameters:*
       - context (Context)
            A Context object used to control the evaluation of integrals
            occurring in the definition of the potential operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the Helmholtz equation
                nabla^2 u + k^2 u = 0.

    *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.

    Note about BEM++ terminology: a *potential operator* acts on functions
    defined on a surface S and produces functions defined at any point of the
    space surrounding S, but not necessarily on S itself. In contrast, a
    *boundary operator* acts on on functions defined on a surface S and produces
    functions defined on the same surface S.
    """
    return _constructHelmholtzPotentialOperator(
        "helmholtz3dFarFieldSingleLayerPotentialOperator", context, waveNumber)

def createHelmholtz3dFarFieldDoubleLayerPotentialOperator(context, waveNumber):
    """
    Create and return a potential operator used to calculate part of the far-field
    pattern for a radiating solution of the Helmholtz equation in 3D.

    See the C++ documentation for mathematical background.

    *Parameters:*
       - context (Context)
            A Context object used to control the evaluation of integrals
            occurring in the definition of the potential operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the Helmholtz equation
                nabla^2 u + k^2 u = 0.

    *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.

    Note about BEM++ terminology: a *potential operator* acts on functions
    defined on a surface S and produces functions defined at any point of the
    space surrounding S, but not necessarily on S itself. In contrast, a
    *boundary operator* acts on on functions defined on a surface S and produces
    functions defined on the same surface S.
    """
    return _constructHelmholtzPotentialOperator(
        "helmholtz3dFarFieldDoubleLayerPotentialOperator", context, waveNumber)

def _constructModifiedHelmholtzOperator(
        className, context,
        domain, range, dualToRange, waveNumber,
        label, useInterpolation, interpPtsPerWavelength):
    basisFunctionType = context.basisFunctionType()
    if (basisFunctionType != domain.basisFunctionType() or
            basisFunctionType != range.basisFunctionType() or
            basisFunctionType != dualToRange.basisFunctionType()):
        raise TypeError("BasisFunctionType of context and all spaces must be the same")
    resultType = context.resultType()

    waveNumberIsComplex = isinstance(waveNumber,complex)
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
    if not label:
        label = ""
    symmetry = 0
    result = _constructObjectTemplatedOnBasisKernelAndResult(
        core, className, basisFunctionType, kernelType, resultType,
        context, domain, range, dualToRange, waveNumber, label,
        symmetry, useInterpolation, interpPtsPerWavelength)
    result._context = context
    result._domain = domain
    result._range = range
    result._dualToRange = dualToRange
    return result

def createModifiedHelmholtz3dSingleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber,
        label=None, useInterpolation=False, interpPtsPerWavelength=5000):
    """
    Create and return a single-layer-potential boundary operator for
    the modified Helmholtz equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the modified Helmholtz equation
                nabla^2 u - k^2 u = 0.
            Only real wave numbers are allowed if context.resultType() is a real
            type (float32 or float64).
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.
       - useInterpolation (bool)
            If set to False (default), the standard exp() function will be used
            to evaluate the exponential factor occurring in the kernel. If set
            to True, the exponential factor will be evaluated by piecewise-cubic
            interpolation of values calculated in advance on a regular
            grid. This normally speeds up calculations, but might result in a
            loss of accuracy. This is an experimental feature: use it at your
            own risk.
       - interpPtsPerWavelength (int)
            If useInterpolation is set to True, this parameter determines the
            number of points per "effective wavelength" (defined as 2 pi /
            abs(waveNumber)) used to construct the interpolation grid. The
            default value (5000) is normally enough to reduce the relative or
            absolute error, *whichever is smaller*, below 100 * machine
            precision. If useInterpolation is set to False, this parameter is
            ignored.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructModifiedHelmholtzOperator(
        "modifiedHelmholtz3dSingleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber,
        label, useInterpolation, interpPtsPerWavelength)

def createModifiedHelmholtz3dDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber,
        label=None, useInterpolation=False, interpPtsPerWavelength=5000):
    """
    Create and return a double-layer-potential boundary operator for
    the modified Helmholtz equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the modified Helmholtz equation
                nabla^2 u - k^2 u = 0.
            Only real wave numbers are allowed if context.resultType() is a real
            type (float32 or float64).
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.
       - useInterpolation (bool)
            If set to False (default), the standard exp() function will be used
            to evaluate the exponential factor occurring in the kernel. If set
            to True, the exponential factor will be evaluated by piecewise-cubic
            interpolation of values calculated in advance on a regular
            grid. This normally speeds up calculations, but might result in a
            loss of accuracy. This is an experimental feature: use it at your
            own risk.
       - interpPtsPerWavelength (int)
            If useInterpolation is set to True, this parameter determines the
            number of points per "effective wavelength" (defined as 2 pi /
            abs(waveNumber)) used to construct the interpolation grid. The
            default value (5000) is normally enough to reduce the relative or
            absolute error, *whichever is smaller*, below 100 * machine
            precision. If useInterpolation is set to False, this parameter is
            ignored.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructModifiedHelmholtzOperator(
        "modifiedHelmholtz3dDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber,
        label, useInterpolation, interpPtsPerWavelength)

def createModifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber,
        label=None, useInterpolation=False, interpPtsPerWavelength=5000):
    """
    Create and return an adjoint double-layer-potential boundary operator for
    the modified Helmholtz equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the modified Helmholtz equation
                nabla^2 u - k^2 u = 0.
            Only real wave numbers are allowed if context.resultType() is a real
            type (float32 or float64).
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.
       - useInterpolation (bool)
            If set to False (default), the standard exp() function will be used
            to evaluate the exponential factor occurring in the kernel. If set
            to True, the exponential factor will be evaluated by piecewise-cubic
            interpolation of values calculated in advance on a regular
            grid. This normally speeds up calculations, but might result in a
            loss of accuracy. This is an experimental feature: use it at your
            own risk.
       - interpPtsPerWavelength (int)
            If useInterpolation is set to True, this parameter determines the
            number of points per "effective wavelength" (defined as 2 pi /
            abs(waveNumber)) used to construct the interpolation grid. The
            default value (5000) is normally enough to reduce the relative or
            absolute error, *whichever is smaller*, below 100 * machine
            precision. If useInterpolation is set to False, this parameter is
            ignored.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructModifiedHelmholtzOperator(
        "modifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber,
        label, useInterpolation, interpPtsPerWavelength)

def createModifiedHelmholtz3dHypersingularBoundaryOperator(
         context, domain, range, dualToRange, waveNumber,
         label=None, useInterpolation=False, interpPtsPerWavelength=5000):
    """
    Create and return a hypersingular boundary operator for the modified
    Helmholtz equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the modified Helmholtz equation
                nabla^2 u - k^2 u = 0.
            Only real wave numbers are allowed if context.resultType() is a real
            type (float32 or float64).
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.
       - useInterpolation (bool)
            If set to False (default), the standard exp() function will be used
            to evaluate the exponential factor occurring in the kernel. If set
            to True, the exponential factor will be evaluated by piecewise-cubic
            interpolation of values calculated in advance on a regular
            grid. This normally speeds up calculations, but might result in a
            loss of accuracy. This is an experimental feature: use it at your
            own risk.
       - interpPtsPerWavelength (int)
            If useInterpolation is set to True, this parameter determines the
            number of points per "effective wavelength" (defined as 2 pi /
            abs(waveNumber)) used to construct the interpolation grid. The
            default value (5000) is normally enough to reduce the relative or
            absolute error, *whichever is smaller*, below 100 * machine
            precision. If useInterpolation is set to False, this parameter is
            ignored.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructModifiedHelmholtzOperator(
        "modifiedHelmholtz3dHypersingularBoundaryOperator", context,
        domain, range, dualToRange, waveNumber,
        label, useInterpolation, interpPtsPerWavelength)

_constructMaxwellOperator = _constructHelmholtzOperator

def _constructModifiedHelmholtzPotentialOperator(className, context, waveNumber):
    basisFunctionType = context.basisFunctionType()
    resultType = context.resultType()
    result = _constructObjectTemplatedOnBasis(
        core, className, basisFunctionType, waveNumber)
    result._context = context
    return result

def createModifiedHelmholtz3dSingleLayerPotentialOperator(context, waveNumber):
    """
    Create and return a single-layer potential operator for the modified Helmholtz
    equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object used to control the evaluation of integrals
            occurring in the definition of the potential operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the modified Helmholtz equation
                -nabla^2 u + k^2 u = 0.

    *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.

    Note about BEM++ terminology: a *potential operator* acts on functions
    defined on a surface S and produces functions defined at any point of the
    space surrounding S, but not necessarily on S itself. In contrast, a
    *boundary operator* acts on on functions defined on a surface S and produces
    functions defined on the same surface S.
    """
    return _constructModifiedHelmholtzPotentialOperator(
        "modifiedHelmholtz3dSingleLayerPotentialOperator", context, waveNumber)

def createModifiedHelmholtz3dDoubleLayerPotentialOperator(context, waveNumber):
    """
    Create and return a double-layer potential operator for the modified Helmholtz
    equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object used to control the evaluation of integrals
            occurring in the definition of the potential operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the modified Helmholtz equation
                -nabla^2 u + k^2 u = 0.

    *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.

    Note about BEM++ terminology: a *potential operator* acts on functions
    defined on a surface S and produces functions defined at any point of the
    space surrounding S, but not necessarily on S itself. In contrast, a
    *boundary operator* acts on on functions defined on a surface S and produces
    functions defined on the same surface S.
    """
    return _constructModifiedHelmholtzPotentialOperator(
        "modifiedHelmholtz3dDoubleLayerPotentialOperator", context, waveNumber)



def createMaxwell3dSingleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber, label=None,
        useInterpolation=False, interpPtsPerWavelength=5000):
    """
    Create and return a single-layer-potential boundary operator for
    the Maxwell equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - waveNumber (float or complex)
            Wave number, i.e. the number k in the Maxwell equation.
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.
       - useInterpolation (bool)
            If set to False (default), the standard exp() function will be used
            to evaluate the exponential factor occurring in the kernel. If set
            to True, the exponential factor will be evaluated by piecewise-cubic
            interpolation of values calculated in advance on a regular
            grid. This normally speeds up calculations, but might result in a
            loss of accuracy. This is an experimental feature: use it at your
            own risk.
       - interpPtsPerWavelength (int)
            If useInterpolation is set to True, this parameter determines the
            number of points per "effective wavelength" (defined as 2 pi /
            abs(waveNumber)) used to construct the interpolation grid. The
            default value (5000) is normally enough to reduce the relative or
            absolute error, *whichever is smaller*, below 100 * machine
            precision. If useInterpolation is set to False, this parameter is
            ignored.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructMaxwellOperator(
        "maxwell3dSingleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber,
        label, useInterpolation, interpPtsPerWavelength)

def createMaxwell3dDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber, label=None,
        useInterpolation=False, interpPtsPerWavelength=5000):
    """
    Create and return a double-layer-potential boundary operator for
    the Maxwell equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - waveNumber (float or complex)
            Wave number.
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.
       - useInterpolation (bool)
            If set to False (default), the standard exp() function will be used
            to evaluate the exponential factor occurring in the kernel. If set
            to True, the exponential factor will be evaluated by piecewise-cubic
            interpolation of values calculated in advance on a regular
            grid. This normally speeds up calculations, but might result in a
            loss of accuracy. This is an experimental feature: use it at your
            own risk.
       - interpPtsPerWavelength (int)
            If useInterpolation is set to True, this parameter determines the
            number of points per "effective wavelength" (defined as 2 pi /
            abs(waveNumber)) used to construct the interpolation grid. The
            default value (5000) is normally enough to reduce the relative or
            absolute error, *whichever is smaller*, below 100 * machine
            precision. If useInterpolation is set to False, this parameter is
            ignored.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructMaxwellOperator(
        "maxwell3dDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber,
        label, useInterpolation, interpPtsPerWavelength)

_constructMaxwellPotentialOperator = _constructHelmholtzPotentialOperator

def createMaxwell3dSingleLayerPotentialOperator(context, waveNumber):
    """
    Create and return a single-layer potential operator for the Maxwell
    equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object used to control the evaluation of integrals
            occurring in the definition of the potential operator.
       - waveNumber (float or complex)
            Wave number.

    *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.

    Note about BEM++ terminology: a *potential operator* acts on functions
    defined on a surface S and produces functions defined at any point of the
    space surrounding S, but not necessarily on S itself. In contrast, a
    *boundary operator* acts on on functions defined on a surface S and produces
    functions defined on the same surface S.
    """
    return _constructMaxwellPotentialOperator(
        "maxwell3dSingleLayerPotentialOperator", context, waveNumber)

def createMaxwell3dDoubleLayerPotentialOperator(context, waveNumber):
    """
    Create and return a double-layer potential operator for the Maxwell
    equation in 3D.

    *Parameters:*
       - context (Context)
            A Context object used to control the evaluation of integrals
            occurring in the definition of the potential operator.
       - waveNumber (float or complex)
            Wave number.

    *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.

    Note about BEM++ terminology: a *potential operator* acts on functions
    defined on a surface S and produces functions defined at any point of the
    space surrounding S, but not necessarily on S itself. In contrast, a
    *boundary operator* acts on on functions defined on a surface S and produces
    functions defined on the same surface S.
    """
    return _constructMaxwellPotentialOperator(
        "maxwell3dDoubleLayerPotentialOperator", context, waveNumber)

def createIdentityOperator(context, domain, range, dualToRange, label=None):
    """
    Create and return a (generalized) identity operator.

    Let X and Y be two function spaces defined on the same grid and represented
    by Space objects supplied in the arguments 'domain' and 'range'. If X is a
    superset of Y (in the mathematical sense), the object returned by this
    function represents the orthogonal projection operator from X to Y. If X is
    a subset of Y, the returned object represents the inclusion operator from X
    to Y. If X is equal to Y, the returned object represents the standard
    identity operator.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.

    All the three spaces ('domain', 'range' and 'dualToRange') must be defined
    on the same grid.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructOperator(
        "identityOperator", context, domain, range, dualToRange, label)

def createLaplaceBeltrami3dOperator(
        context, domain, range, dualToRange, label=None):
    """
    Create and return a Laplace-Beltrami operator in 3D.

    The weak form of this operator is

        int_S (grad_S phi^* . grad_S psi) dS

    where grad_S is the surface gradient operator, ^* denotes complex
    conjugation, and phi and psi are the test and trial functions, respectively.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.

    All the three spaces ('domain', 'range' and 'dualToRange') must be defined
    on the same grid.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructOperator(
        "laplaceBeltrami3dOperator", context, domain, range, dualToRange, label)

def createMaxwell3dFarFieldSingleLayerPotentialOperator(context, waveNumber):
    """
    Create and return a potential operator used to calculate part of the
    far-field pattern for a radiating solution of the Maxwell equations in 3D.

    See the C++ documentation for mathematical background.

    *Parameters:*
       - context (Context)
            A Context object used to control the evaluation of integrals
            occurring in the definition of the potential operator.
       - waveNumber (float or complex)
            Wave number.

    *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.

    Note about BEM++ terminology: a *potential operator* acts on functions
    defined on a surface S and produces functions defined at any point of the
    space surrounding S, but not necessarily on S itself. In contrast, a
    *boundary operator* acts on on functions defined on a surface S and produces
    functions defined on the same surface S.
    """
    return _constructMaxwellPotentialOperator(
        "maxwell3dFarFieldSingleLayerPotentialOperator", context, waveNumber)

def createMaxwell3dFarFieldDoubleLayerPotentialOperator(context, waveNumber):
    """
    Create and return a potential operator used to calculate part of the
    far-field pattern for a radiating solution of the Maxwell equations in 3D.

    See the C++ documentation for mathematical background.

    *Parameters:*
       - context (Context)
            A Context object used to control the evaluation of integrals
            occurring in the definition of the potential operator.
       - waveNumber (float or complex)
            Wave number.

    *Returns* a newly constructed PotentialOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.

    Note about BEM++ terminology: a *potential operator* acts on functions
    defined on a surface S and produces functions defined at any point of the
    space surrounding S, but not necessarily on S itself. In contrast, a
    *boundary operator* acts on on functions defined on a surface S and produces
    functions defined on the same surface S.
    """
    return _constructMaxwellPotentialOperator(
        "maxwell3dFarFieldDoubleLayerPotentialOperator", context, waveNumber)

def createMaxwell3dIdentityOperator(context, domain, range, dualToRange,
                                    label=None):
    """
    Create and return an identity operator with weak form defined under an
    antisymmetric pseudo-inner product.

    This class represents an operator I_M whose weak form (under the standard
    inner product) is

    <u, I_M v> = \int_S u^* . (v \times n) dS,

    where n is the outward unit vector normal to the surface S at a given point.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.

    All the three spaces ('domain', 'range' and 'dualToRange') must be defined
    on the same grid.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructOperator(
        "maxwell3dIdentityOperator", context, domain, range, dualToRange, label)

def createNullOperator(context, domain, range, dualToRange, label=None):
    """
    Create and return a null (zero-valued) operator.

    *Parameters:*
       - context (Context)
            A Context object to control the assembly of the weak form of the
            newly constructed operator.
       - domain (Space)
            Function space to be taken as the domain of the operator.
       - range (Space)
            Function space to be taken as the range of the operator.
       - dualToRange (Space)
            Function space to be taken as the dual to the range of the operator.
       - label (string)
            Textual label of the operator. If set to None (default), a unique
            label will be generated automatically.

    *Returns* a newly constructed BoundaryOperator_BasisFunctionType_ResultType
    object, with BasisFunctionType and ResultType determined automatically from
    the context argument and equal to either float32, float64, complex64 or
    complex128.
    """
    return _constructOperator(
        "nullOperator", context, domain, range, dualToRange, label)

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
        core, "Python" + functorType,
        resultType, function, argumentDimension, resultDimension)
    result = _constructObjectTemplatedOnBasisAndResult(
        core, "gridFunctionFromPython" + functorType,
        basisFunctionType, resultType,
        context, space, dualSpace, functor)
    result._context = context
    result._space = space
    result._dualSpace = dualSpace
    return result

def createGridFunction(
        context, space, dualSpace=None,
        function=None, surfaceNormalDependent=False, coefficients=None, projections=None):
    """
    Create and return a GridFunction object with values determined by a Python
    function or by an input vector of coefficients or projections.

    *Parameters:*
       - context (Context)
            Assembly context from which a quadrature strategy can be retrieved.
       - space (Space)
            Function space to expand the grid function in.
       - dualSpace (Space)
            Function space dual to 'space'.
       - coefficients (Vector)
            A vector of the coefficients of the function in the basis of space
            'space'.
       - projections (Vector)
            A vector of projections of the function on the basis of space
            'dualSpace'.
       - function (a Python callable object)
            Function object whose values on 'space.grid()' will be used to
            construct the new grid function. If 'surfaceNormalDependent' is set
            to False (default), 'function' will be passed a single argument
            containing a 1D array of the coordinates of a point lying on the
            grid 'space.grid()'. If 'surfaceNormalDependent' is set to True, the
            function will be passed one more argument, a 1D array containing the
            components of the unit vector normal to 'space.grid()' at the point
            given in the first argument. In both cases 'function' should return
            its value at the given point, in the form of a scalar or a 1D array
            with dimension equal to 'space.codomainDimension()'.
       - surfaceNormalDependent (bool)
            Indicates whether the grid function depends on the unit vector
            normal to the grid or not.

    If both the 'space' and 'dualSpace' are given, they must be defined on the
    same grid and have the same codomain dimension.

    This function can be called with three different sets of parameters.

    Variant 1: construction of a grid function from a Python function::

        createGridFunction(context, space, dualSpace, function,
                           surfaceNormalDependent)

    Variant 2: construction of a grid function from the vector of its
    coefficients in the basis of the space 'space'::

        createGridFunction(context, space, coefficients=coefficients)

    Variant 3: construction of a grid function from the vector of its
    projections on the basis functions of the space 'dualSpace'::

        createGridFunction(context, space, dualSpace, projections=projections)

    Example scalar-valued function defined in a 3D space that can be passed to
    'createGridFunction' with 'surfaceNormalDependent = False'::

        def fun1(point):
            x, y, z = point
            r = np.sqrt(x**2 + y**2 + z**2)
            return 2 * x * z / r**5 - y / r**3

    Example scalar-valued function defined in a 3D space that can be passed to
    'createGridFunction' with 'surfaceNormalDependent = True'::

        import cmath
        def fun2(point, normal):
            x, y, z = point
            nx, ny, nz = normal
            k = 5
            return cmath.exp(1j * k * x) * (nx - 1)
    """

    params = [function,coefficients,projections]
    params_active = sum([0 if p is None else 1 for p in params])
    if params_active != 1 :
        raise ValueError("createGridFunction(): Exactly one of 'function', "
                         "'coefficients' or 'projections' must be supplied")
    if function is not None and dualSpace is None:
        raise ValueError("createGridFunction(): You must set the dualSpace "
                         "parameter to a valid Space object when constructing "
                         "a grid function from a Python function")
    if projections is not None and dualSpace is None:
        raise ValueError("createGridFunction(): You must set the dualSpace "
                         "parameter to a valid Space object when constructing "
                         "a grid function from a vector of projections")

    if coefficients is not None:
        return _constructObjectTemplatedOnBasisAndResult(
            core, "gridFunctionFromCoefficients", context.basisFunctionType(),
            context.resultType(), context, space, dualSpace, coefficients)
    if projections is not None:
        return _constructObjectTemplatedOnBasisAndResult(
            core, "gridFunctionFromProjections", context.basisFunctionType(),
            context.resultType(), context, space, dualSpace, projections)

    if surfaceNormalDependent:
        className = "SurfaceNormalDependentFunctor"
    else:
        className = "SurfaceNormalIndependentFunctor"
    return __gridFunctionFromFunctor(
        className, context, space, dualSpace, function,
        argumentDimension=space.grid().dimWorld(),
        resultDimension=space.codomainDimension())

def gridFunctionFromSurfaceNormalDependentFunction(
        context, space, dualSpace, function,
        argumentDimension=3, resultDimension=1):
    """
    Deprecated. Superseded by createGridFunction().
    """
    print ("gridFunctionFromSurfaceNormalDependentFunction(): DEPRECATED. "
           "Please use the createGridFunction() function instead.")
    return __gridFunctionFromFunctor(
        "SurfaceNormalDependentFunctor",
        context, space, dualSpace, function,
        argumentDimension, resultDimension)

def gridFunctionFromSurfaceNormalIndependentFunction(
        context, space, dualSpace, function,
        argumentDimension=3, resultDimension=1):
    """
    Deprecated. Superseded by createGridFunction().
    """
    print ("gridFunctionFromSurfaceNormalDependentFunction(): DEPRECATED. "
           "Please use the createGridFunction() function instead.")
    return __gridFunctionFromFunctor(
        "SurfaceNormalIndependentFunctor",
        context, space, dualSpace, function,
        argumentDimension, resultDimension)

def estimateL2Error(
        gridFunction, exactFunction, quadStrategy, evaluationOptions,
        surfaceNormalDependent=False):
    """
    Calculate the L^2-norm of the difference between a grid function and a function
    defined as a Python callable.

    *Parameters:*
       - gridFunction (GridFunction)
            A grid function.
       - exactFunction (a Python callable object)
            A Python callable object whose values on 'gridFunction.grid()' will
            be compared against those of the grid function. If
            'surfaceNormalDependent' is set to False (default), 'exactFunction'
            will be passed a single argument containing a 1D array of the
            coordinates of a point lying on the grid 'gridFunction.grid()'. If
            'surfaceNormalDependent' is set to True, the function will be passed
            one more argument, a 1D array containing the components of the unit
            vector normal to 'gridFunction.grid()' at the point given in the
            first argument. In both cases 'exactFunction' should return its
            value at the given point, in the form of a scalar or a 1D array with
            dimension equal to 'space.codomainDimension()'.
       - quadStrategy (QuadratureStrategy)
            Quadrature stategy to be used to evaluate integrals involved in the
            L^2-norm.
       - evaluationOptions (EvaluationOptions)
            Additional options controlling function evaluation.
       - surfaceNormalDependent (bool)
            Indicates whether 'exactFunction' depends on the unit vector
            normal to the grid or not.

    *Returns* a tuple (absError, relError), where absError is the L^2-norm of
    the difference between 'gridFunction' and 'exactFunction' and relError is
    equal to absError divided by the L^2-norm of 'exactFunction'

    See the documentation of createGridFunction() for example definitions of
    Python functions that can be passed to the 'exactFunction' argument.
    """
    basisFunctionType = checkType(gridFunction.basisFunctionType())
    resultType = checkType(gridFunction.resultType())
    if (basisFunctionType != quadStrategy.basisFunctionType() or
            resultType != quadStrategy.resultType()):
        raise TypeError("BasisFunctionType and ResultType of gridFunction "
                        "and exactFunction must be the same")
    if surfaceNormalDependent:
        dependent = "Dependent"
    else:
        dependent = "Independent"
    functor = _constructObjectTemplatedOnValue(
        core, "PythonSurfaceNormal%sFunctor" % dependent,
        resultType, exactFunction,
        gridFunction.space().grid().dimWorld(), # argument dimension
        gridFunction.space().codomainDimension() # result dimension
        )
    absError = _constructObjectTemplatedOnBasisAndResult(
        core, "L2NormOfDifferenceFromPythonSurfaceNormal%sFunctor" % dependent,
        basisFunctionType, resultType,
        gridFunction, functor, quadStrategy, evaluationOptions)
    exactFunctionL2Norm = _constructObjectTemplatedOnBasisAndResult(
        core, "L2NormOfDifferenceFromPythonSurfaceNormal%sFunctor" % dependent,
        basisFunctionType, resultType,
        0 * gridFunction, functor, quadStrategy, evaluationOptions)
    return absError, absError / exactFunctionL2Norm

def createDefaultIterativeSolver(
        boundaryOperator,
        convergenceTestMode="test_convergence_in_dual_to_range"):
    """
    Create and return a DefaultIterativeSolver object.

    The DefaultIterativeSolver class acts as an interface to Belos, the
    iterative solver package from Trilinos. It lets you solve the equation A f =
    g for the function f, with A being a boundary operator and g a grid
    function.

    *Parameters:*
       - boundaryOperator (BoundaryOperator or BlockedBoundaryOperator)
            The boundary operator A standing on the left-hand-side of the
            equation to be solved.
       - convergenceTestMode (string)
            Convergence test mode. Can be either
            "test_convergence_in_dual_to_range" (default) or
            "test_convergence_in_range". See below.

    Convergence can be tested either in the range space of the operator A or in
    the space dual to the range. A standard Galerkin discretisation of the form
    Ax = b maps into the space dual to the range of the operator. If you choose
    to test in the range space, the equation pinv(M)Ax = pinv(M)b is solved,
    where M is the mass matrix mapping from the range space into its dual and
    pinv(M) is its pseudoinverse.

    *Returns* a newly constructed
    DefaultIterativeSolver_BasisFunctionType_ResultType object, with
    BasisFunctionType and ResultType determined automatically from the
    boundaryOperator argument and equal to either float32, float64, complex64 or
    complex128.
    """
    basisFunctionType = boundaryOperator.basisFunctionType()
    resultType = boundaryOperator.resultType()
    result = _constructObjectTemplatedOnBasisAndResult(
        core, "DefaultIterativeSolver", basisFunctionType, resultType,
        boundaryOperator, convergenceTestMode)
    result._boundaryOperator = boundaryOperator
    return result

def createDefaultDirectSolver(
        boundaryOperator,
        convergenceTestMode="test_convergence_in_dual_to_range"):
    """
    Create and return a DefaultDirectSolver object.

    The DefaultDirectSolver class acts as an interface to the
    direct solver from LAPACK. It lets you solve the equation A f =
    g for the function f, with A being a boundary operator and g a grid
    function.

    *Parameters:*
       - boundaryOperator (BoundaryOperator or BlockedBoundaryOperator)
            The boundary operator A standing on the left-hand-side of the
            equation to be solved.

    *Returns* a newly constructed
    DefaultDirectSolver_BasisFunctionType_ResultType object, with
    BasisFunctionType and ResultType determined automatically from the
    boundaryOperator argument and equal to either float32, float64, complex64 or
    complex128.
    """
    basisFunctionType = boundaryOperator.basisFunctionType()
    resultType = boundaryOperator.resultType()
    result = _constructObjectTemplatedOnBasisAndResult(
        core, "DefaultDirectSolver", basisFunctionType, resultType,
        boundaryOperator)
    result._boundaryOperator = boundaryOperator
    return result

from core import defaultGmresParameterList
from core import defaultCgParameterList

def createAccuracyOptions():
    """Create and return an AccuracyOptions object with default settings."""
    return core.AccuracyOptions()

def createAccuracyOptionsEx():
    """Create and return an AccuracyOptionsEx object with default settings."""
    return core.AccuracyOptionsEx()

def createAcaOptions():
    """Create and return an AcaOptions object with default settings."""
    return core.AcaOptions()

def createEvaluationOptions():
    """Create and return an EvaluationOptions object with default settings."""
    return core.EvaluationOptions()

def createBlockedOperatorStructure(context):
    """
    Create and return a BlockedOperatorStructure object.

    BlockedOperatorStructure is a helper class used in construction of blocked
    boundary operators. It represents a matrix of boundary operators. To
    construct a blocked boundary operator, store individual boundary operators
    in appropriate rows and columns of a BlockedOperatorStructure object,
    repeatedly calling its setBlock() method, and pass this object to the
    createBlockedBoundaryOperator() function for validation.

    *Parameters:*
       - context (Context)
            A Context object. The values returned by
            context.basisFunctionType() and context.resultType() will
            determine the precise type of the newly constructed
            BlockedOperatorStructure object.

    *Returns* a newly constructed
    BlockedOperatorStructure_BasisFunctionType_ResultType object, with
    BasisFunctionType and ResultType determined automatically from the 'context'
    argument and equal to either float32, float64, complex64 or complex128.
    """
    name = 'BlockedOperatorStructure'
    return _constructObjectTemplatedOnBasisAndResult(
        core, name, context.basisFunctionType(), context.resultType())

def createBlockOperatorStructure(context):
    """
    Deprecated. Superseded by createBlockedOperatorStructure().
    """
    print ("createBlockOperatorStructure(): DEPRECATED. Please use the "
           "createBlockedOperatorStructure() function instead.")
    return createBlockedOperatorStructure(context)

def createBlockedBoundaryOperator(context, structure):
    """
    Create and return a BlockedBoundaryOperator object.

    A blocked boundary operator is an operator that consists of several blocks
    arranged in a matrix, each of which is an "elementary" boundary operator
    represented by a BoundaryOperator object.

    *Parameters:*
       - structure (BlockedOperatorStructure or nested list)

            An object determining the boundary operators to be put in specific
            blocks of the newly constructed blocked boundary operator. Can be
            either a BlockedBoundaryOperator object or simply a nested (2D)
            list. In the latter case, empty blocks can be marked with None.

    All the boundary operators from a single column of 'structure' must have the
    same domain, and all the operators from a single row of 'structure' must have
    the same range and space dual to range. No row and no column of 'structure'
    must be completely empty (contain only uninitialized BoundaryOperator
    objects). If these conditions are not satisfied, an exception is thrown.

    *Returns* a newly constructed
    BlockedBoundaryOperator_BasisFunctionType_ResultType object, with
    BasisFunctionType and ResultType determined automatically from the
    'structure' argument and equal to either float32, float64, complex64 or
    complex128.

    *Example:*

    Assume that A, B and C are boundary operators. We want to create
    the blocked operator ::

            [A B]
        Z = [C 0],

    where 0 is an empty block.

    Variant 1::

        from bempp import lib
        context = lib.createContext(...)
        structure = lib.createBlockedOperatorStructure(context)
        structure.setBlock(0, 0, A)
        structure.setBlock(0, 1, B)
        structure.setBlock(1, 0, C)
        Z = lib.createBlockedBoundaryOperator(context, structure)

    Variant 2::

        from bempp import lib
        context = lib.createContext(...)
        Z = lib.createBlockedBoundaryOperator(context, [[A, B], [C, None]])
    """
    structureClass = getattr(
        core, "BlockedOperatorStructure_" + context.basisFunctionType() +
        "_" + context.resultType())
    if isinstance(structure, structureClass):
        print "real structure"
        boStructure = structure
    else:
        boStructure = createBlockedOperatorStructure(context)
        for i in range(len(structure)):
            for j in range(len(structure[i])):
                block = structure[i][j]
                if block:
                    boStructure.setBlock(i, j, block)
    name = 'BlockedBoundaryOperator'
    return _constructObjectTemplatedOnBasisAndResult(
        core, name, context.basisFunctionType(), context.resultType(), boStructure)

def adjoint(*args):
    """
    Return the adjoint of a BoundaryOperator.

    *Note* The BoundaryOperator must act on real-valued basis functions,
    otherwise an exception is thrown.
    """
    return _constructObjectTemplatedOnBasisAndResult(
        core, "adjoint", args[0].basisFunctionType(), args[0].resultType(), *args)

if core._withAhmed:
    def acaOperatorApproximateLuInverse(operator, delta):
        """
        Create and return a discrete boundary operator representing an approximate
        inverse of an H-matrix.

        *Parameters:*
           - operator (DiscreteBoundaryOperator)
                A discrete boundary operator stored in the form of an H-matrix.
           - delta (float)
                Approximation accuracy.

        *Returns* a DiscreteBoundaryOperator_ValueType object representing an
        approximate inverse of the operator supplied in the 'operator' argument,
        stored in the form of an approximate H-matrix LU decomposition. ValueType is
        set to operator.valueType().
        """
        # name = 'createAcaApproximateLuInverse'
        name = 'acaOperatorApproximateLuInverse'
        return _constructObjectTemplatedOnValue(
            core, name, operator.valueType(), operator, delta)

    def createAcaApproximateLuInverse(operator, delta):
        """
        Deprecated. Superseded by acaOperatorApproximateLuInverse().
        """
        print ("createAcaApproximateLuInverse(): DEPRECATED: "
               "use acaOperatorApproximateLuInverse() instead.")
        return acaOperatorApproximateLuInverse(operator, delta)

    def scaledAcaOperator(operator, multiplier):
        """
        Multiply a discrete boundary operator stored as an H-matrix by a scalar.

        *Parameters:*
           - operator (DiscreteBoundaryOperator)
                A discrete boundary operator stored in the form of an H-matrix.
           - multiplier (float or complex, depending on operator.valueType())
                Scalar with which the supplied operator should be multiplied.

        *Returns* a newly constructed DiscreteBoundaryOperator_ValueType object
        storing an H-matrix equal to the H-matrix stored in 'operator' and
        multiplied by 'multiplier'. ValueType is set to operator.valueType().
        """
        name = 'scaledAcaOperator'
        return _constructObjectTemplatedOnValue(
            core, name, operator.valueType(), operator, multiplier)

    def acaOperatorSum(op1, op2, eps, maximumRank):
        """
        Create and return a discrete boundary operator representing an approximate
        sum of two discrete boundary operators stored as H-matrices.

        *Parameters:*
           - op1 (DiscreteBoundaryOperator)
                First operand; a discrete boundary operator stored in the form of an
                H-matrix.
           - op2 (DiscreteBoundaryOperator)
                Second operand; a discrete boundary operator stored in the form of an
                H-matrix.
           - eps (float)
                Approximation accuracy. (TODO: explain better)
           - maximumRank (int)
                Maximum rank of blocks that should be considered low-rank in the
                H-matrix to be constructed.

        *Returns* a newly constructed DiscreteBoundaryOperator_ValueType object
        storing an H-matrix approximately equal to sum of the H-matrices stored in
        the two operands. ValueType is set to op1.valueType() (which must be equal
        to op2.valueType()).
        """
        name = 'acaOperatorSum'
        if (op1.valueType() != op2.valueType()):
            raise TypeError("acaOperatorSum: ValueTypes of 'op1' and 'op2' do not match.")
        return _constructObjectTemplatedOnValue(
            core, name, op1.valueType(), op1, op2, eps, maximumRank)

def discreteDenseBoundaryOperator(mat):
    """
    Return a discrete operator object stored as a dense matrix mat.

    *Parameters:*
       - mat
           A two-dimensional array to be wrapped by the discrete operator.

    *Returns* a newly constructed DiscreteBoundaryOperator_ValueType object
    wrapping the matrix mat.
    """
    name = 'discreteDenseBoundaryOperator'
    mat = numpy.asanyarray(mat)
    dtype = numpy.promote_types(numpy.float32, mat.dtype).name
    return _constructObjectTemplatedOnValue(core, name, dtype, mat)

def discreteSparseInverse(op):
    """
    Return a discrete operator object that evaluates the inverse of op applied
    to a vector by constructing the LU decomposition. For this operation op must
    represent a sparse matrix operator.

    *Parameters:*
       - op (DiscreteBoundaryOperator)
           A discrete boundary operator stored as sparse matrix

    *Returns* a newly constructed DiscreteBoundaryOperator_ValueType object
    which stores the LU decomposition of op and evaluates the inverse of op
    applied to a vector.
    """
    name = 'discreteSparseInverse'
    return _constructObjectTemplatedOnValue(core, name, op.valueType(), op)

def discreteOperatorToPreconditioner(op):
    """
    Create a preconditioner from a discrete boundary operator.

    *Parameters:*
       - op (DiscreteBoundaryOperator)
           A discrete operator, which acts as preconditioner, e.g. an ACA
           approximate LU decomposition or a sparse inverse.

    *Returns* a preconditioner.
    """
    name = 'discreteOperatorToPreconditioner'
    return _constructObjectTemplatedOnValue(core, name, op.valueType(), op)

def discreteBlockDiagonalPreconditioner(opArray):
    """
    Create a block diagonal preconditioner from a list of operators.

    *Parameters:*
       - opArray (list of DiscreteBoundaryOperator objects)
           A list of DiscreteBoundaryOperator objects op1, op2, ..., opn.

    *Returns* a block diagonal preconditioner P = diag(op1, op2, ..., opn).
    """
    name = 'discreteBlockDiagonalPreconditioner'
    return _constructObjectTemplatedOnValue(
        core, name, opArray[0].valueType(), opArray)

def areInside(grid, points):
    """Determine whether points lie inside a grid (assumed to be closed)."""
    if points.shape[0] != grid.dimWorld():
        raise ValueError("insideClosedGrid(): the number of columns "
                         "in the array 'points' must match the number of "
                         "dimensions of the world in which 'grid' is embedded")
    return numpy.array(core.areInside(grid, points))

from bempp.core import GridSegment, gridSegmentWithPositiveX
