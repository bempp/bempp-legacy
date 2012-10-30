#!/usr/bin/python

# This script postprocesses the HTML documentation for the Python
# interface of BEM++. In particular, it makes the docs generic
# with respect to BasisFunctionType and ResultType.

import errno, glob, os, re, shutil, sys

html_dir = sys.argv[1]
if not os.path.isdir(html_dir):
    raise IOError(errno.ENOENT, "'" + html_dir + "' is not a directory")
paths = glob.glob(os.path.join(html_dir, "*.html"))

print "Postprocessing Python documentation..."

entities_templated_on_BFT = [
    "Space",
    "piecewiseConstantScalarSpace"
    "piecewiseLinearContinuousScalarSpace"
    ]
entities_templated_on_VT = [
    "DiscreteBoundaryOperator",
    "Preconditioner",
    "PythonSurfaceNormalIndependentFunctor",
    "PythonSurfaceNormalDependentFunctor",
    "acaBlockDiagonalPreconditioner",
    "acaOperatorApproximateLuInverse",
    "acaOperatorSum",
    "scaledAcaOperator",
    "discreteSparseInverse",
    "enable_shared_from_this_discrete_boundary_operator",
    "LinearOpDefaultBase"
    ]
entities_templated_on_RT = [
    "helmholtz3dAdjointDoubleLayerBoundaryOperator",
    "helmholtz3dDoubleLayerBoundaryOperator",
    "helmholtz3dHypersingularBoundaryOperator",
    "helmholtz3dSingleLayerBoundaryOperator",
    "helmholtz3dDoubleLayerPotentialOperator",
    "helmholtz3dSingleLayerPotentialOperator",
    "helmholtz3dFarFieldDoubleLayerPotentialOperator",
    "helmholtz3dFarFieldSingleLayerPotentialOperator",
    ]
entities_templated_on_BFT_and_RT = [
    "AbstractBoundaryOperator",
    "BlockedBoundaryOperator",
    "BlockedOperatorStructure",
    "BlockedSolution",
    "BoundaryOperator",
    "Context",
    "DefaultIterativeSolver",
    "GridFunction",
    "PotentialOperator",
    "QuadratureStrategyBase",
    "QuadratureStrategy",
    "SolutionBase",
    "Solution",
    "Solver",
    "identityOperator",
    "nullOperator",
    "laplace3dAdjointDoubleLayerBoundaryOperator",
    "laplace3dDoubleLayerBoundaryOperator",
    "laplace3dHypersingularBoundaryOperator",
    "laplace3dSingleLayerBoundaryOperator",
    "numericalQuadratureStrategy",
    "gridFunctionFromCoefficients",
    "gridFunctionFromPythonSurfaceNormalDependentFunctor",
    "gridFunctionFromPythonSurfaceNormalIndependentFunctor"
    ]
entities_templated_on_BFT_KT_and_RT = [
    "modifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator",
    "modifiedHelmholtz3dDoubleLayerBoundaryOperator",
    "modifiedHelmholtz3dHypersingularBoundaryOperator",
    "modifiedHelmholtz3dSingleLayerBoundaryOperator",
]

re_python_BFT = re.compile(
    r'(?<!id="bempp\.core\.)(?<!#bempp\.core\.)' +
    "(" + "|".join(entities_templated_on_BFT) + ")" +
    "_complex128"
    )
replacement_python_BFT = r"\1_<i>BasisFunctionType</i>"
re_cpp_BFT = re.compile(
    r"(?:Fiber|Bempp)::" +
    "(" + "|".join(entities_templated_on_BFT) + ")" +
    r"&lt;"
    r"\(std::complex&lt;\(double\)&gt;\)"
    r"&gt;"
    )
replacement_cpp_BFT = r"\1&lt;<i>BasisFunctionType</i>&gt;"

re_python_RT = re.compile(
    r'(?<!id="bempp\.core\.)(?<!#bempp\.core\.)' +
    "(" + "|".join(entities_templated_on_RT) + ")" +
    "_complex128"
    )
replacement_python_RT = r"\1_<i>ResultType</i>"
re_cpp_RT = re.compile(
    r"(?:Fiber|Bempp)::" +
    "(" + "|".join(entities_templated_on_RT) + ")" +
    r"&lt;"
    r"\(std::complex&lt;\(double\)&gt;\)"
    r"&gt;"
    )
replacement_cpp_RT = r"\1&lt;<i>ResultType</i>&gt;"

re_python_VT = re.compile(
    r'(?<!id="bempp\.core\.)(?<!#bempp\.core\.)' +
    "(" + "|".join(entities_templated_on_VT) + ")" +
    "_complex128"
    )
replacement_python_VT = r"\1_<i>ValueType</i>"
re_cpp_VT = re.compile(
    r"(?:Fiber|Bempp)::" +
    "(" + "|".join(entities_templated_on_VT) + ")" +
    r"&lt;"
    r"\(std::complex&lt;\(double\)&gt;\)"
    r"&gt;"
    )
replacement_cpp_VT = r"\1&lt;<i>ValueType</i>&gt;"


re_python_BFT_and_RT = re.compile(
    r'(?<!id="bempp\.core\.)(?<!#bempp\.core\.)(?<!DONT_REPLACE)' +
    "(" + "|".join(entities_templated_on_BFT_and_RT) + ")" +
    "_complex128_complex128"
    )
replacement_python_BFT_and_RT = r"\1_<i>BasisFunctionType</i>_<i>ResultType</i>"
re_cpp_BFT_and_RT = re.compile(
    r"(?:Fiber|Bempp)::" +
    "(" + "|".join(entities_templated_on_BFT_and_RT) + ")" +
    r"&lt;"
    r"\(std::complex&lt;\(double\)&gt;"
    ","
    r"std::complex&lt;\(double\)&gt;\)"
    r"&gt;"
    )
replacement_cpp_BFT_and_RT = r"\1&lt;<i>BasisFunctionType</i>,<i>ResultType</i>&gt;"


re_python_BFT_KT_and_RT = re.compile(
    r'(?<!id="bempp\.core\.)(?<!#bempp\.core\.)' +
    "(" + "|".join(entities_templated_on_BFT_KT_and_RT) + ")" +
    "_complex128_complex128_complex128"
    )
replacement_python_BFT_KT_and_RT = r"\1_<i>BasisFunctionType</i>_<i>KernelType</i>_<i>ResultType</i>"
re_cpp_BFT_KT_and_RT = re.compile(
    r"(?:Fiber|Bempp)::" +
    "(" + "|".join(entities_templated_on_BFT_KT_and_RT) + ")" +
    r"&lt;"
    r"\(std::complex&lt;\(double\)&gt;"
    ","
    r"std::complex&lt;\(double\)&gt;"
    ","
    r"std::complex&lt;\(double\)&gt;\)"
    r"&gt;"
    )
replacement_cpp_BFT_KT_and_RT = r"\1&lt;<i>BasisFunctionType</i>,<i>KernelType</i>,<i>ResultType</i>&gt;"

for path in paths:
    f = open(path, "r")
    text = f.read()
    f.close()

    text = re_python_BFT_KT_and_RT.sub(replacement_python_BFT_KT_and_RT, text)
    text = re_cpp_BFT_KT_and_RT.sub(replacement_cpp_BFT_KT_and_RT, text)

    text = re_python_BFT_and_RT.sub(replacement_python_BFT_and_RT, text)
    text = re_cpp_BFT_and_RT.sub(replacement_cpp_BFT_and_RT, text)

    text = re_python_BFT.sub(replacement_python_BFT, text)
    text = re_cpp_BFT.sub(replacement_cpp_BFT, text)

    text = re_python_RT.sub(replacement_python_RT, text)
    text = re_cpp_RT.sub(replacement_cpp_RT, text)

    text = re_python_VT.sub(replacement_python_VT, text)
    text = re_cpp_VT.sub(replacement_cpp_VT, text)

    text = text.replace("Fiber::QuadratureStrategy&lt;(std::complex&lt;(double)&gt;,"
                        "std::complex&lt;(double)&gt;,Bempp::GeometryFactory)&gt;",
                        "Fiber::QuadratureStrategy&lt;<i>BasisFunctionType</i>,"
                        "<i>ResultType</i>,Bempp::GeometryFactory)&gt;")
    text = text.replace("Fiber::QuadratureStrategyBase&lt;(std::complex&lt;(double)&gt;,"
                        "std::complex&lt;(double)&gt;,Bempp::GeometryFactory)&gt;",
                        "Fiber::QuadratureStrategyBase&lt;<i>BasisFunctionType</i>,"
                        "<i>ResultType</i>,Bempp::GeometryFactory)&gt;")
    text = text.replace("DONT_REPLACE", "")

    shutil.copy2(path, path + ".orig")
    f = open(path, "w")
    f.write(text)
    f.close()


