Steklov-Poincare Operator
=========================

In this tutorial we're going to use BEM++ to implement an approximation to a Dirichlet to Neumann 
operator for the exterior Laplace problem. 

Mathematical background
-----------------------
Suppose that a (smooth) domain, :math:`\Omega \subset \mathbb R^3` has boundary :math:`\Gamma = \partial \Omega` and 
that :math:`u` is a function satisfying :math:`\triangle u = 0` in :math:`\Omega^c`.  The exterior Dirichlet and Neumann traces, 
:math:`\gamma_0^{ext} u` and :math:`\gamma_1^{ext} u` satisfy

.. math:: V \gamma_1^{ext} u = (-\frac{1}{2} I + K) \gamma_0^{ext} u
    :label: steklov

where :math:`V` and :math:`K` are the single and double layer operators associated with the fundamental solution of 
the Laplace operator.    

We suppose that :math:`\{\phi_i, i=1\dots n\}` is a basis for a function space on :math:`\Gamma` and that we 
know the Dirichlet data :math:`\gamma_0^{ext}u =: g = \sum_i g_i \phi_i`.  We will seek an approximation
to the Neumann data, :math:`\gamma_1^{ext}u \approx b = \sum_i b_i \phi_i` such that :eq:`steklov` is satisfied weakly, i.e.

.. math:: \sum_i g_i (V \phi_i, \phi_j) = \sum_i b_i((-\frac{1}{2} I + K)\phi_i, \phi_j) \qquad \forall j 

Implementation
--------------


.. blahmath::
    V:H^{-\frac{1}{2}}(\Gamma) \rightarrow H^{\frac{1}{2}}(\Gamma)
    
    K:H^{\frac{1}{2}}(\Gamma) \rightarrow H^{\frac{1}{2}}(\Gamma)