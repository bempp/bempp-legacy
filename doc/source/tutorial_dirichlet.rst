Steklov-Poincare Operator
=========================

In this tutorial we're going to use BEM++ to implement an approximation to a Dirichlet to Neumann 
operator for the exterior Laplace problem. 

Mathematical background
-----------------------
Suppose that a (smooth) domain, :math:`\Omega \subset \mathbb R^3` has boundary :math:`\Gamma = \partial \Omega` and 
that :math:`u` is a function satisfying :math:`\triangle u = 0` in :math:`\Omega^c` that decays to zero at infinity.  
The exterior Dirichlet and Neumann traces, :math:`\gamma_0^{ext} u` and :math:`\gamma_1^{ext} u` satisfy

.. math:: V \gamma_1^{ext} u = (-\frac{1}{2} I + K) \gamma_0^{ext} u
    :label: steklov

where 

.. math::
    V:H^{-\frac{1}{2}}(\Gamma) \rightarrow H^{\frac{1}{2}}(\Gamma)
    
    K:H^{\frac{1}{2}}(\Gamma) \rightarrow H^{\frac{1}{2}}(\Gamma)

are the single and double layer operators associated with the fundamental solution of 
the Laplace operator and :math:`I` is the identity.      

Suppose that we have function spaces, :math:`S^0 \subset H^{-\frac{1}{2}}(\Gamma)` and 
:math:`S^1 \subset H^{\frac{1}{2}}(\Gamma)` with bases :math:`\{\phi^0_i, i=1\dots n_0\}` 
and :math:`\{\phi^1_i, i=1\dots n_1\}` respectively and that we know the 
know the Dirichlet data 

.. math:: \gamma_0^{ext}u =: g = \sum_i g_i \phi^0_i  

We will seek an approximation to the Neumann data, 

.. math:: \gamma_1^{ext}u \approx b = \sum_i b_i \phi^1_i 

such that :eq:`steklov` is satisfied weakly, i.e.

.. math:: \sum_i g_i (V \phi^0_i, \phi^0_j) = \sum_i b_i((-\frac{1}{2} I + K)\phi^1_i, \phi^0_j) \qquad \forall j 

Note that the :math:`H^{-\frac{1}{2}}(\Gamma)`-ellipticity of :math:`V` means that the problem is well-posed. 
In this tutorial, we will use BEM++ to assemble the matrices, :math:`V_{ij} := (V \phi^0_i, \phi^0_j)` and 
:math:`M_{ij} := ((-\frac{1}{2} I + K)\phi^1_i, \phi^0_j)`.  

.. highlight:: c

Implementation
--------------
A complete listing can be found in ``examples/tutorial_dirichlet.cpp``.  

First we need to define the domain, these are described by a ``Bempp::Grid`` object in BEM++. To import a GMSH grid, run::

   #include "grid/grid_factory.hpp"
   using namespace Bempp;
 
   ...
   
   std::string myGmshFilename = "sphere-644.msh";
   GridParameters params;
   params.topology = GridParameters::TRIANGULAR;
   std::auto_ptr<Grid> grid = GridFactory::importGmshGrid(params, myGmshFilename, false, false);
   
Now we can define some approximation spaces.  For this example, we will use a piecewise linear space for :math:`S^1` 
and piecewice constants for :math:`S^0`::

   #include "space/piecewise_linear_continuous_scalar_space.hpp"
   #include "space/piecewise_constant_scalar_space.hpp"
   
   ...
   
   PiecewiseLinearContinuousScalarSpace<double> S0(*grid);
   PiecewiseConstantScalarSpace<double> S1(*grid);
   
   S1.assignDofs();
   S0.assignDofs();

The calls to assignDofs() initialise the spaces.

Next we need to build some operators.  BEM++ has the concept of two types of operator.  The underlying "continuous"
linear operators are represented by lightweight classes that implement the ``Bempp::LinearOperator`` interface.

::

   #include "assembly/identity_operator.hpp"
   #include "assembly/single_layer_potential_3d.hpp"
   #include "assembly/double_layer_potential_3d.hpp"

   ...
   
   SingleLayerPotential3D<double> slp;
   DoubleLayerPotential3D<double> dlp;
   IdentityOperator<double> id;

To perform calculations with these operators, we need to discretise them.  Discretised operators are subclasses
of ``Bempp::DiscreteScalarValuedLinearOperator``.  To perform the discretisation, we need to evaluate (some of) the terms
:math:`V_{ij}` and :math:`M_{ij}`.  This is done using Fiber.   

::

   #include "fiber/standard_local_assembler_factory_for_operators_on_surfaces.hpp"
   #include "fiber/standard_local_assembler_factory_for_operators_on_surfaces.hpp"
   
   ...
   
   Fiber::AccuracyOptions accuracyOptions; // default
   Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory> factory(accuracyOptions);

The Fiber assembly factory provides the LinearOperators with assembly objects which are able to perform 
the integrations.  All we need to do is pass it to them:
   
::   

   #include "assembly/assembly_options.hpp"
   #include "assembly/discrete_scalar_valued_linear_operator.hpp"
   
   ...
   
   AssemblyOptions assemblyOptions;
   typedef std::auto_ptr<DiscreteScalarValuedLinearOperator<double> > DiscreteLinearOperatorPtr;
   DiscreteLinearOperatorPtr discreteSlp =
      slp.assembleWeakForm(pwConstSpace, pwLinearCtsSpace, factory, assemblyOptions);
   DiscreteLinearOperatorPtr discreteDlp =
      dlp.assembleWeakForm(pwConstSpace, pwLinearCtsSpace, factory, assemblyOptions);
   DiscreteLinearOperatorPtr discreteId =
      id.assembleWeakForm(pwConstSpace, pwLinearCtsSpace, factory, assemblyOptions);

The ``Bempp::AssemblyOptions`` and ``Fiber::AccuracyOptions`` and classes allow us to configure the type of
discretisation and integration.  In this case, we are content with the default options.  BEM++ supports
the representation of the discrete operator in several forms.  In this example we're just going to ask
for dense local Matrices, which are managed using the Armadillo C++ linear algebra library. 

::

   #include <armadillo>

   ...
   
   arma::Mat<double> M = -0.5 * discreteId->asMatrix() + discreteDlp->asMatrix();
   arma::Mat<double> V = discreteSlp->asMatrix();
 
Finally, we can solve the system. Both the bases that we have used are nodal, so determining the :math:`g_i` 
and interpreting the :math:`b_i` is straightforward.  In this case, we're imposing constant Dirichlet data.  
On the unit sphere contained in sphere-644.msh, this means that the exterior Laplace solution is :math:`u(x) = \frac{1}{|x|}`,
so we expect the Neumann data to be uniformly equal to -1 

::

   arma::Col<double> g = arma::ones(linearCtsSpace.globalDofCount(), 1);
   arma::Col<double> b = arma::solve(V, M * g);
   std::cout<<"Neumann coefficients"<<b;
