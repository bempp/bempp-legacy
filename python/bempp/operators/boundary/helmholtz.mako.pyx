#cython: embedsignature=True

__all__=['single_layer','double_layer','adjoint_double_layer','hypersingular', 'osrc_dtn']

from bempp.utils.parameter_list cimport ParameterList
from bempp.space.space cimport Space
from bempp.operators.boundary import modified_helmholtz as _modified_helmholtz
from bempp.assembly import BoundaryOperatorBase
import numpy as np

<% ops = [('single_layer','Return the Helmholtz single layer boundary operator.'),
          ('double_layer','Return the Helmholtz double layer boundary operator.'),
          ('adjoint_double_layer','Return the Helmholtz adjoint double layer boundary operator.'),
          ('hypersingular','Return the Helmholtz hypersingular boundary operator.')]
    
%>

% for op,help_text in ops:
def ${op}(Space domain, Space range, Space dual_to_range,
        object wave_number,
        object label="",
        object symmetry="auto_symmetry",
        object parameters=None):
    """ 

    ${help_text}

    """

    return _modified_helmholtz.${op}(domain,range,dual_to_range,
            wave_number/1j,label,symmetry,parameters)


% endfor

def osrc_dtn(space, wave_number, npade = 2, theta = np.pi/3):
    return _OsrcDtn(space, wave_number, npade, theta)

class _OsrcDtn(BoundaryOperatorBase):

    def __init__(self, space, wave_number, npade = 2, theta = np.pi/3.):

        # cdef object _space
        # cdef object _wave_number
        # cdef object _basis_type
        # cdef object _result_type
        # cdef object _npade

        self._space = space
        self._wave_number = wave_number
        self._result_type = 'complex128'
        self._basis_type = 'float64'
        self._npade = npade
        self._theta = theta

        self._weak_form = None

    def weak_form(self):

        from bempp.operators.boundary.sparse import identity
        from bempp.operators.boundary.sparse import laplace_beltrami
        from bempp import ZeroBoundaryOperator
        from bempp import InverseSparseDiscreteBoundaryOperator

        if self._weak_form is not None:
            return self._weak_form

        space = self._space

        # Get the operators
        mass = identity(space, space, space).weak_form()
        stiff = laplace_beltrami(space, space, space).weak_form()

        # Compute damped wavenumber

        bbox = self._space.grid.bounding_box
        rad = np.linalg.norm(bbox[1,:]-bbox[0,:])/2
        dk = self._wave_number+ 1j * 0.4 * self._wave_number**(1.0/3.0) * rad**(-2.0/3.0)

        # Get the Pade coefficients

        c0, alpha, beta = _pade_coeffs(self._npade, self._theta)

        # Now put all operators together

        term = ZeroBoundaryOperator(space, space, space).weak_form()

        for i in range(self._npade):
            term = term - (alpha[i]/(dk**2) * stiff * 
                    InverseSparseDiscreteBoundaryOperator(
                        mass - beta[i]/(dk**2) * stiff))
        result = 1j*self._wave_number*(c0*mass + term * mass)

        self._weak_form = result
        
        return result
        
    @property
    def domain(self):
        return self._space

    @property
    def range(self):
        return self._space

    @property
    def dual_to_range(self):
        return self._space


    
def osrc_ntd(space, wave_number, npade = 2, theta = np.pi/3):
    return _OsrcNtd(space, wave_number, npade, theta)

class _OsrcNtd(BoundaryOperatorBase):

    def __init__(self, space, wave_number, npade = 2, theta = np.pi/3.):

        # cdef object _space
        # cdef object _wave_number
        # cdef object _basis_type
        # cdef object _result_type
        # cdef object _npade

        self._space = space
        self._wave_number = wave_number
        self._result_type = 'complex128'
        self._basis_type = 'float64'
        self._npade = npade
        self._theta = theta

        self._weak_form = None

    def weak_form(self):

        from bempp.operators.boundary.sparse import identity
        from bempp.operators.boundary.sparse import laplace_beltrami
        from bempp import ZeroBoundaryOperator
        from bempp import InverseSparseDiscreteBoundaryOperator

        if self._weak_form is not None:
            return self._weak_form

        space = self._space

        # Get the operators
        mass = identity(space, space, space).weak_form()
        stiff = laplace_beltrami(space, space, space).weak_form()

        # Compute damped wavenumber

        bbox = self._space.grid.bounding_box
        rad = np.linalg.norm(bbox[1,:]-bbox[0,:])/2
        dk = self._wave_number+ 1j * 0.4 * self._wave_number**(1.0/3.0) * rad**(-2.0/3.0)

        # Get the Pade coefficients

        c0, alpha, beta = _pade_coeffs(self._npade, self._theta)

        # Now put all operators together

        term = ZeroBoundaryOperator(space, space, space).weak_form()

        for i in range(self._npade):
            term = term - (alpha[i]/(dk**2) * stiff * 
                    InverseSparseDiscreteBoundaryOperator(
                        mass - beta[i]/(dk**2) * stiff))
        result = 1./(1j*self._wave_number)*(
                mass * InverseSparseDiscreteBoundaryOperator(mass-1./(dk**2)*stiff)*(c0*mass + term * mass))

        self._weak_form = result
        
        return result
        
    @property
    def domain(self):
        return self._space

    @property
    def range(self):
        return self._space

    @property
    def dual_to_range(self):
        return self._space

def _pade_coeffs(n, theta):
    """
    _pade_coeffs(n, theta)
    
    Compute the coefficients of the Pade approximation of (1 + Delta_Gamma/k_eps^2)^0.5,
    for a given order of expansion and the angle of the branch cut.
    See X. Antoine et al., CMAME 195 (2006).
    See M. Darbas et al., J. Comp. Physics 236 (2013).
    """


    # First define the coefficients for the standard Pade approximation.
    Np = n

    #C0 = 1.0
    aj = np.zeros(Np)
    bj = np.zeros(Np)
    for jj in range(1,Np+1): # remember Python skips the last index
        # remember that Python starts with zero in an array
        aj[jj-1] = 2.0 / (2.0 * Np + 1.0) * np.sin(jj * np.pi / (2.0 * Np + 1.0))**2
        bj[jj-1] = np.cos(jj * np.pi / (2.0 * Np + 1.0))**2
    # Now define the coefficients for the 'theta' branch cut.
    C0t = np.exp(1j * theta / 2.0) * (1.0 + np.sum( ( aj * (np.exp(-1j * theta) - 1.0 ) ) / ( 1.0 + bj * (np.exp(-1j * theta) - 1.0 ) ) ) )
    ajt = np.exp(-1j * theta / 2.0) * aj / ( (1 + bj * (np.exp(-1j*theta) - 1))**2 )
    bjt = np.exp(-1j * theta) * bj / (1 + bj * (np.exp(-1j*theta) - 1))
    return C0t, ajt, bjt





        


