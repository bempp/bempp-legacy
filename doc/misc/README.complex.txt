Date: 27.04.2012. Updated: 6.09.2012
Author: Wojtek

On the support for complex-valued operators and basis functions.

There are four new CMake options:

ENABLE_SINGLE_PRECISION,
ENABLE_DOUBLE_PRECISION,
ENABLE_COMPLEX_KERNELS,
ENABLE_COMPLEX_BASIS_FUNCTIONS;

I hope their names are self-explanatory. By default (currently) all are ON. They
determine which explicit template instantiations take place during compilation.

Specific template classes depend on different sets of parameters. In particular:

* The BoundaryOperator template class depends on two parameters: BasisFunctionType
(type representing individual components of the functions making up the test and
trial spaces used to discretise the operator) and ResultType (type of the
elements of the discretised operator's matrix).

* The GridFunction template class also depends on BasisFunctionType (type used
to represent individual components of the basis functions from the space in
which the grid function is expanded) and ResultType (type of the expansion
coefficients).

* The DiscreteLinearOperator depends only on the type of its matrix entries,
denoted ValueType.

* The Solver, a little non-intuitively, depends on BasisFunctionType and
ResultType. This is because it acts on LinearOperators and GridFunctions, which
depend on these parameters.

* Several classes from the Fiber module (like the integrators) depend in
addition on KernelType (type used to represent the values of a kernel function).
However, they are instantiated only for those combinations of BasisFunctionType,
KernelType and ResultType that make sense. Users are not expected to need to
access these classes.

In addition, many classes internally typedef a CoordinateType, which is
guaranteed to be a real type.

Some issues for further consideration:

1. At present, objects with different ResultType cannot be used together (e.g.
arithmetic operations between BoundaryOperators with different ResultType are not
defined; similarly, a linear operator cannot act on a GridFunction with a
different ResultType). This means in particular that integral operators with
complex-valued kernels (and hence a complex ResultType) need to be used in
combination with identity operators with complex ResultType. I don't think this
is a serious problem.

2. The explicit dependence of high-level classes such as BoundaryOperator and
GridFunction on BasisFunctionType is somewhat unsatisfactory, since this type
reveals itself in their public interface only very weakly (via the
Space<BasisFunctionType> type of arguments passed into their constructors). We
could remove this dependence by overloading the constructors to accept spaces of
both real- and complex-valued functions. This would make the interface nicer,
but it would also strongly increase the implementation complexity (basically the
problem is that C++ doesn't support templated virtual functions, hence we would
need to manually overload all virtual functions taking an object parametrised by
BasisFunctionType, such as a Space or a Basis). In my opinion, it isn't worth
the effort.

3. Template parameter names could perhaps be chosen more judiciously. In
particular, ResultType is a bit too vague. BasisFunctionValueType would probably
also be more clear than BasisFunctionType, but it is rather long and unwieldy.
Ideas are welcome!
