////////////////////////////////////////////////////////////////////////////////
// File: classBempp_1_1AbstractBoundaryOperator.xml
%define AbstractBoundaryOperator_docstring(BASIS, RESULT) 
"Abstract (non-discretized) boundary operator.

An AbstractBoundaryOperator represents a linear mapping $L : X \\\\to
Y$ between two function spaces $X : S \\\\to K^p$ (_domain_) and $Y :
T \\\\to Q^q$ (_range_) defined on $n$-dimensional surfaces $S$ and
$T$ embedded in an $(n+1)$-dimensional domain. Each of the symbols
$K$ and $Q$ can stand either for the set of real or complex
numbers. The surfaces $S$ and $T$ may be equal.

The function assembleWeakForm() can be used to construct the weak form
of the operator.

The values of the (components of the) basis functions into which the
functions from $X$ are expanded are represented with objects of type 
B"
#BASIS ". The values of the (components of the) functions from $Y$,
i.e. the results of operator's action, are represented with objects of
type " #RESULT "."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    AbstractBoundaryOperator)

// /*  Construction and destruction  */

// %feature("docstring")
// Bempp::AbstractBoundaryOperator::AbstractBoundaryOperator "

// Constructor.

// Parameters:
// -----------

// domain:  Function space being the domain of the operator.

// range:  Function space being the range of the operator.

// dualToRange:  Function space dual to the the range of the operator.

// label:  Textual label of the operator. If empty, a unique label is
// generated automatically.

// symmetry:  Symmetry of the weak form of the operator. Can be any
// combination of the flags defined in the enumeration type Symmetry.

// None of the shared pointers may be null and the spaces range and
// dualToRange must be defined on the same grid, otherwise an exception
// is thrown. ";

/*  Spaces  */

%define AbstractBoundaryOperator_domain_docstring(BASIS, RESULT) 
"Return the function space being the domain of this operator."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    AbstractBoundaryOperator, domain)

%define AbstractBoundaryOperator_range_docstring(BASIS, RESULT) 
"Return the function space being the range of this operator."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    AbstractBoundaryOperator, range)

%define AbstractBoundaryOperator_dualToRange_docstring(BASIS, RESULT) 
"Return the function space dual to the range of this operator."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    AbstractBoundaryOperator, dualToRange)

/*  Other attributes  */

%define AbstractBoundaryOperator_label_docstring(BASIS, RESULT) 
"Return the label of the operator."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    AbstractBoundaryOperator, label)

%define AbstractBoundaryOperator_symmetry_docstring(BASIS, RESULT) 
"Return the symmetry properties of the operator.

The returned value should be treated as a bitwise combination of the
flags defined in the Symmetry enumeration type.

TODO: convert this into a human-readable string."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    AbstractBoundaryOperator, symmetry)

%define AbstractBoundaryOperator_isLocal_docstring(BASIS, RESULT) 
"Return whether this operator is local.

Suppose that an operator $A$ acting on a function $f(x)$ produces
another function $g(x)$. We say that $A$ is local if the value of $g$
at any point $x$ depends only on the values of $f$ in an infinitesimal
neighbourhood of $x$.

Multiplicative and differential operators are local and discretization
of their weak forms with finite elements leads to sparse matrices.
Conversely, integral operators are in general non-local and
discretization of their weak forms leads to dense matrices. "
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    AbstractBoundaryOperator, isLocal)

/*  Assembly  */

%define AbstractBoundaryOperator_assembleWeakForm_docstring(BASIS, RESULT) 
"Assemble and return the operator's weak form.

This function constructs a discrete linear operator representing the
matrix $L_{jk}$ with entries of the form

\\\\[L_{jk} = \\\\int_S \\\\phi_j L \\\\psi_k,\\\\]

where $L$ is the linear operator represented by this object, $S$
denotes the surface on which the domain function space $X$ is defined
and which is represented by the grid returned by domain.grid(),
$\\\\phi_j$ is a _test function_ from the space $Y'$ dual to the range
of the operator, $Y$, and $\\\\psi_k$ is a _trial function_ from the
domain space $X$. "
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    AbstractBoundaryOperator, assembleWeakForm)


////////////////////////////////////////////////////////////////////////////////
// File: classBempp_1_1BoundaryOperator.xml
%define BoundaryOperator_docstring(BASIS, RESULT) 
"Operator acting on functions defined on a surface.

A BoundaryOperator is a lightweight wrapper of an
AbstractBoundaryOperator and a DiscreteBoundaryOperator representing
its weak form. The weak form is evaluated lazily, on the
first call to weakForm(). The Context object passed to the constructor
of the BoundaryOperator or to the initialize() function determines how
this weak form is calculated.

Different threads should not share BoundaryOperator objects, since
notably the weakForm() function is not thread-safe. Instead, each
thread should hold its own copy of a BoundaryOperator (note that
copying BoundaryOperators is cheap -- the copy constructor is
shallow). See the documentation of AbstractBoundaryOperator for
information on the meaning of the types B" #BASIS " and R" #RESULT "."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    BoundaryOperator)

// %define BoundaryOperator_BoundaryOperator_docstring(BASIS)
// "Construct an uninitialized BoundaryOperator. ";

// %define BoundaryOperator_BoundaryOperator_docstring(BASIS)

// Construct and initialize a BoundaryOperator.

// Equivalent to calling the initialize() function on a BoundaryOperator
// object created with the default constructor. See the documentation of
// initialize() for a description of the constructor's parameters.

// User code should not need to invoke this constructor directly; it is
// more convenient to use non-member constructors supplied with
// particular AbstractBoundaryOperator subclasses (e.g.
// laplace3dSingleLayerBoundaryOperator(), identityOperator(),
// pseudoinverse() etc.), which construct an AbstractBoundaryOperator and
// wrap it in a BoundaryOperator in a single step. ";

%define BoundaryOperator_initialize_docstring(BASIS, RESULT)
"Initialize or reinitialize a BoundaryOperator.

*Arguments*
    context (Context_" #BASIS "_" #RESULT ")
        Context that will be used to build the weak form of abstractOp
        when necessary.
    abstractOp (AbstractBoundaryOperator_" #BASIS "_" #RESULT ")
        AbstractBoundaryOperator that will
        be encapsulated in this BoundaryOperator.

The provided objects are stored internally. In addition, any
previously calculated weak form of the abstract boundary operator is
invalidated."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    BoundaryOperator, initialize)

%define BoundaryOperator_uninitialize_docstring(BASIS, RESULT)
"Uninitialize a BoundaryOperator.

This function erases any internal references to a context, an abstract
boundary operator and its weak form."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    BoundaryOperator, uninitialize)

%define BoundaryOperator_isInitialized_docstring(BASIS, RESULT)
"Return true if the BoundaryOperator has been initialized, false otherwise."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    BoundaryOperator, isInitialized)

%define BoundaryOperator_abstractOperator_docstring(BASIS, RESULT)
"Return the encapsulated abstract boundary operator
or None if the operator has not been initialized yet."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    BoundaryOperator, abstractOperator)

%define BoundaryOperator_context_docstring(BASIS, RESULT)
"Return the Context object with which the operator has been initialized,
or None if the operator has not been initialized yet."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    BoundaryOperator, context)

%define BoundaryOperator_weakForm_docstring(BASIS, RESULT)
"Return the weak form of the encapsulated abstract boundary operator.

An exception is thrown if this function is called on an uninitialized
BoundaryOperator."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    BoundaryOperator, weakForm)

%define BoundaryOperator_domain_docstring(BASIS, RESULT)
"Return the domain of the encapsulated abstract boundary operator, or
None if the operator has not been initialized yet."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    BoundaryOperator, domain)

%define BoundaryOperator_range_docstring(BASIS, RESULT)
"Return the range of the encapsulated abstract boundary operator, or
None if the operator has not been initialized yet."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    BoundaryOperator, range)

%define BoundaryOperator_dualToRange_docstring(BASIS, RESULT)
"Return the space dual to the range of the encapsulated abstract boundary
operator, or None if the operator has not been initialized yet."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    BoundaryOperator, dualToRange)

%define BoundaryOperator_label_docstring(BASIS, RESULT)
"Return the label of this BoundaryOperator."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    BoundaryOperator, label)

%define BoundaryOperator_apply_docstring(BASIS, RESULT)
"Act on a grid function.

*Arguments*

This function sets y_inout to alpha * A * x_in + beta * y_inout, where
A is the operator represented by this object.

The space of x_in must be identical with the domain of the
encapsulated abstract boundary operator, whereas the space of y_inout
and its dual must be identical with the range of the encapsulated
abstract boundary operator and its dual; otherwise an exception is
thrown. An exception is also thrown if the BoundaryOperator is
uninitialized."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    BoundaryOperator, apply)


/////////////////////////////////////////////////////////////////////////////////
// File: classBempp_1_1Context.xml
%define Context_docstring(BASIS, RESULT) 
"Assembly context.

This class manages the assembly of weak forms and evaluation of
potentials.

An assembly context consists of a quadrature strategy, which
determines the way integrals are calculated, and assembly options,
which control higher-level aspects of weak-form assembly, e.g. the use
or not of acceleration algorithms such as ACA and the level of
parallelism.

A call to the Context::getWeakForm() function returns the weak form of
the abstract boundary operator passed in the argument to this
function, assembled using the settings from the given Context. The
Context may store the weak form in an internal cache; see the
documentation of getWeakForm() for more details.

See the documentation of AbstractBoundaryOperator for information
about the meaning of the types B" #BASIS " and R" #RESULT "."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    Context)

// %define Context_Context_docstring(BASIS, RESULT)

// Constructor.

// Parameters:
// -----------

// quadStrategy:  Quadrature strategy to be used for calculation of
// integrals occurring e.g. in the weak forms of boundary operators or in
// the definition of potential operators.

// assemblyOptions:  Further options influencing the weak-form assembly
// process. ";

%define Context_getWeakForm_docstring(BASIS, RESULT) 
"Return the weak form of the specified abstract operator.

This function returns returns the weak form of the specified abstract
boundary operator, calculated in accordance with the settings
specified during the construction of the Context.

An important design principle in BEM++ is that abstract boundary
operators are immutable after construction. The same is true for
Context objects. Therefore, a Context stores newly calculated weak
forms in an internal cache, unless a given abstract boundary operator
does not provide a valid unique identifier. As long as the weak form
remains in cache, subsequent calls to Context::getWeakForm() with the
same abstract boundary operator will not recalculate the weak form,
but return the cached instance. It should, however, be noted that the
cache does not maintain a persistent relationship with weak forms: it
stores them as weak rather than shared pointers. Therefore, a weak
form is deallocated as soon as the last reference to it *apart from
that residing in the cache* goes out of scope.

Note: This function may be removed from the Python interface of BEM++
in a future version of the library."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    Context, getWeakForm)

%define Context_assemblyOptions_docstring(BASIS, RESULT) 
"Return the AssemblyOptions object passed when
constructing the Context. "
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    Context, assemblyOptions)

%define Context_quadStrategy_docstring(BASIS, RESULT) 
"Return the QuadratureStrategy object passed when
constructing the Context. "
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    Context, quadStrategy)

////////////////////////////////////////////////////////////////////////////////
// File: classBempp_1_1Space.xml
%define Space_docstring(BASIS)
"Function space.

This class represents a space of functions defined on a grid. The
space is spanned by a finite number of scalar- or vector-valued basis
functions. The values of (the components of) these basis functions are
represented with objects of type B" #BASIS ".

The basis functions of a space, also known as global degrees of
freedom (DOFs), can have support extending over multiple elements of
the grid. They are, however, composed of one or more *local* basis
functions (local degrees of freedom), each of which resides on a
single element. The mapping of local to global degrees of freedom is
triggered by calling the function assignDofs(). Many other member
functions of Space may only be invoked after assignDofs() has beeen
called."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_CLASS_TEMPLATED_ON_BASIS(
    Space)

/*  Attributes  */

%define Space_domainDimension_docstring(BASIS)
"Return the dimension of the grid on which functions from this space 
are defined."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS(
    Space, domainDimension)

%define Space_codomainDimension_docstring(BASIS)
"Return the dimension of the codomain of the functions.

In other words, number of components of the values of the functions.
(E.g. H1 space -> 1, H(curl) space on a 2D surface -> 2)."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS(
    Space, codomainDimension)

%define Space_grid_docstring(BASIS)
"Return the grid on which the functions from this space are defined."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS(
    Space, grid)

/*  DOF management  */

%define Space_assignDofs_docstring(BASIS)
"Assign global degrees of freedom to local degrees of freedom."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS(
    Space, assignDofs)

%define Space_dofsAssigned_docstring(BASIS)
"Return True if assignDofs() has been called before, False otherwise."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS(
    Space, dofsAssigned)

%define Space_flatLocalDofCount_docstring(BASIS)
"Return the total number of local degrees of freedom on all elements."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS(
    Space, flatLocalDofCount)

%define Space_globalDofCount_docstring(BASIS)
"Return number of global degrees of freedom or zero
if assignDofs() has not been called before."
%enddef
BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS(
    Space, globalDofCount)

// %define Space_Space_docstring(BASIS)

// Constructor.

// Parameters:
// -----------

// grid:   Grid on which functions from this space should be defined.

// The supplied grid must remain valid until the Space object is
// destructed.

// Todo The grid should be passed via a shared pointer instead of via a
// reference. ";

