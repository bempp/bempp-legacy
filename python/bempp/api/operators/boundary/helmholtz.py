# pylint: disable-msg=too-many-arguments
"""Definition of the Helmholtz boundary operators."""
from bempp.api.assembly.boundary_operator import \
    BoundaryOperator as _BoundaryOperator
import numpy as _np


def single_layer(domain, range_, dual_to_range,
                 wave_number,
                 label="SLP", symmetry='no_symmetry',
                 parameters=None,
                 use_projection_spaces=True,
                 assemble_only_singular_part=False):
    """Return the Helmholtz single-layer boundary operator.

    Parameters
    ----------
    domain : bempp.api.space.Space
        Domain space.
    range_ : bempp.api.space.Space
        Range space.
    dual_to_range : bempp.api.space.Space
        Dual space to the range space.
    wave_number : complex
        Wavenumber for the Helmholtz problem.
    label : string
        Label for the operator.
    symmetry : string
        Symmetry mode. Possible values are: 'no_symmetry',
        'symmetric', 'hermitian'.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given the
        default global parameter object `bempp.api.global_parameters`
        is used.
    use_projection_spaces : bool
        Represent operator by projection from higher dimensional space
        if available. This parameter can speed up fast assembly routines,
        such as H-Matrices or FMM (default true).
    assemble_only_singular_part : bool
        When assembled the operator will only contain components for adjacent or 
        overlapping test and trial functions (default false).

    """
    from .modified_helmholtz import single_layer as sl

    return sl(domain, range_, dual_to_range,
              wave_number / (1j), label, symmetry,
              parameters, use_projection_spaces, assemble_only_singular_part)


def double_layer(domain, range_, dual_to_range,
                 wave_number,
                 label="DLP", symmetry='no_symmetry',
                 parameters=None,
                 use_projection_spaces=True,
                 assemble_only_singular_part=False):
    """Return the Helmholtz double-layer boundary operator.

    Parameters
    ----------
    domain : bempp.api.space.Space
        Domain space.
    range_ : bempp.api.space.Space
        Range space.
    dual_to_range : bempp.api.space.Space
        Dual space to the range space.
    wave_number : complex
        Wavenumber for the Helmholtz problem.
    label : string
        Label for the operator.
    symmetry : string
        Symmetry mode. Possible values are: 'no_symmetry',
        'symmetric', 'hermitian'.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given the
        default global parameter object `bempp.api.global_parameters`
        is used.
    use_projection_spaces : bool
        Represent operator by projection from higher dimensional space
        if available. This parameter can speed up fast assembly routines,
        such as H-Matrices or FMM (default true).
    assemble_only_singular_part : bool
        When assembled the operator will only contain components for adjacent or 
        overlapping test and trial functions (default false).

    """

    from .modified_helmholtz import double_layer as dl

    return dl(domain, range_, dual_to_range,
              wave_number / (1j), label, symmetry,
              parameters,
              use_projection_spaces,
              assemble_only_singular_part)


def adjoint_double_layer(domain, range_, dual_to_range,
                         wave_number,
                         label="ADJ_DLP", symmetry='no_symmetry',
                         parameters=None,
                         use_projection_spaces=True,
                         assemble_only_singular_part=False):
    """Return the Helmholtz adjoint double-layer boundary operator.

    Parameters
    ----------
    domain : bempp.api.space.Space
        Domain space.
    range_ : bempp.api.space.Space
        Range space.
    dual_to_range : bempp.api.space.Space
        Dual space to the range space.
    wave_number : complex
        Wavenumber for the Helmholtz problem.
    label : string
        Label for the operator.
    symmetry : string
        Symmetry mode. Possible values are: 'no_symmetry',
        'symmetric', 'hermitian'.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given the
        default global parameter object `bempp.api.global_parameters`
        is used.
    use_projection_spaces : bool
        Represent operator by projection from higher dimensional space
        if available. This parameter can speed up fast assembly routines,
        such as H-Matrices or FMM (default true).
    assemble_only_singular_part : bool
        When assembled the operator will only contain components for adjacent or 
        overlapping test and trial functions (default false).

    """

    from .modified_helmholtz import adjoint_double_layer as adl

    return adl(domain, range_, dual_to_range,
               wave_number / (1j), label, symmetry,
               parameters,
               use_projection_spaces=True,
               assemble_only_singular_part=False)


def hypersingular(domain, range_, dual_to_range,
                  wave_number,
                  label="HYP", symmetry='no_symmetry',
                  parameters=None,
                  use_slp=False,
                  use_projection_spaces=True,
                  assemble_only_singular_part=False):
    """Return the Helmholtz hypersingular boundary operator.

    Parameters
    ----------
    domain : bempp.api.space.Space
        Domain space.
    range_ : bempp.api.space.Space
        Range space.
    dual_to_range : bempp.api.space.Space
        Dual space to the range space.
    wave_number : complex
        Wavenumber for the Helmholtz problem.
    label : string
        Label for the operator.
    symmetry : string
        Symmetry mode. Possible values are: 'no_symmetry',
        'symmetric', 'hermitian'.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given the
        default global parameter object `bempp.api.global_parameters`
        is used.
    use_projection_spaces : bool
        Represent operator by projection from higher dimensional space
        if available. This parameter can speed up fast assembly routines,
        such as H-Matrices or FMM (default true).
    use_slp : True/False or boundary operator object
        The hypersingular operator can be represented as a sparse transformation
        of a single-layer operator. If `use_slp=True` this representation is used.
        If `use_slp=op` for a single-layer boundary operator assembled on a
        suitable space this operator is used to assemble the hypersingular operator.
        Note that if `use_slp=op` is used no checks are performed if the slp operator
        is correctly defined for representing the hypersingular operator. Hence,
        if no care is taken this option can lead to a wrong operator. Also,
        `use_slp=True` or `use_slp=op` is only valid if the `domain` and `dual_to_range`
        spaces are identical.
    assemble_only_singular_part : bool
        When assembled the operator will only contain components for adjacent or 
        overlapping test and trial functions (default false).
        Note. This option is only used if `use_slp` is not specified.
    """

    from .modified_helmholtz import hypersingular as hyp

    return hyp(domain, range_, dual_to_range,
               wave_number / (1j), label, symmetry,
               use_slp=use_slp, parameters=parameters,
               use_projection_spaces=use_projection_spaces,
               assemble_only_singular_part=False)


def multitrace_operator(grid, wave_number, parameters=None, spaces='linear', target=None):
    """Return the Helmholtz multitrace operator.

    Parameters
    ----------
    grid : bempp.api.grid.Grid
        The underlying grid for the multitrace operator
    wave_number : complex
        Wavenumber of the operator.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.
    spaces: string
        Choose 'linear' to assemble the operator
        with continuous linear function spaces for the
        Dirichlet and Neumann component (default). For
        a dual pairing of a linear space for the Dirichlet
        data and piecewise constant space for the Neumann
        data choose 'dual'.
    target: bempp.api.grid.Grid
        Specifies a target grid. If it is different from
        'grid' then the operator maps from 'grid' to
        'target'.

    """

    import bempp.api
    return bempp.api.operators.boundary.modified_helmholtz.multitrace_operator(
        grid, wave_number / (1j), parameters, spaces, target=target)


def single_layer_and_hypersingular_pair(grid, wave_number, parameters=None, spaces='linear', base_slp=None,
                                        return_base_slp=False):
    """Return a pair of hypersingular and single layer operator.

    This function creates a pair of a single-layer and a hypersingular
    operator, where both operators are instantiated using a common
    base single-layer operator. Hence, only one single-layer operator
    needs to be discretized to obtain both operators on the given
    grid.

    Parameters
    ----------
    grid : bempp.api.grid.Grid
        The underlying grid for the multitrace operator
    wave_number : complex
        Wavenumber of the operator.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.
    spaces: string
        Choose 'linear' to assemble the operator
        with continuous linear function spaces for the
        Dirichlet and Neumann component (default). For
        a dual pairing of a linear space for the Dirichlet
        data and piecewise constant space for the Neumann
        data choose 'dual'.
    base_slp : None
        Specify a base single-layer operator to be used. If
        set to None, a base single-layer operator will be
        instantiated by the function.
    return_base_slp : bool
        If True also return the original large space single layer
        operator from which the hypersingular and slp operator
        are derived. Default is False

    Returns
    -------
    A pair (slp, hyp) of a single-layer and hypersingular operator.
    If return_base_slp is true a triplet (slp, hyp, base_slp) is
    returned, where base_slp is the single-layer operator, from
    which slp and hyp are obtained via sparse transformations.

    """
    import bempp.api
    return bempp.api.operators.boundary.modified_helmholtz.single_layer_and_hypersingular_pair(
        grid, wave_number / (1j), parameters, spaces, base_slp, return_base_slp)


def interior_calderon_projector(grid, wave_number, parameters=None):
    """Return the interior Calderon projector.

    Parameters
    ----------
    grid : bempp.api.grid.Grid
        The underlying grid for the multitrace operator
    wave_number : complex
        Wavenumber of the operator.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.

    """

    from .sparse import multitrace_identity

    return .5 * multitrace_identity(grid, parameters) + multitrace_operator(grid, wave_number, parameters)


def exterior_calderon_projector(grid, wave_number, parameters=None):
    """Return the exterior Calderon projector.

    Parameters
    ----------
    grid : bempp.api.grid.Grid
        The underlying grid for the multitrace operator
    wave_number : complex
        Wavenumber of the operator.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.

    """

    from .sparse import multitrace_identity

    return .5 * multitrace_identity(grid, parameters) - multitrace_operator(grid, wave_number, parameters)


def osrc_dtn(space, wave_number, npade=2, theta=_np.pi / 3, damped_wavenumber=None,
             parameters=None, label="osrc_dtn"):
    """Return the OSRC approximation to the DtN map.

    Parameters
    ----------
    space : bempp.api.space.Space
        Underlying function space. Needs to be a continuous space.
    wave_number : complex
        Wavenumber of the operator.
    npade : int
        Number of Pade terms to be used in the approximation
    theta : float64
        Angle of the branch cut of the square root operator.
    damped_wavenumber : complex
        Damped wavenumber of the OSRC operator.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.
    label : string
        Label for the operator.

    """

    return _OsrcDtn(space, wave_number, npade, theta, damped_wavenumber, label=label)


class _OsrcDtn(_BoundaryOperator):

    def __init__(self, space, wave_number, npade=2, theta=_np.pi / 3., damped_wavenumber=None,
                 parameters=None, label="osrc_dtn"):
        super(_OsrcDtn, self).__init__(space, space, space, label=label)

        self._space = space
        self._wave_number = wave_number
        self._npade = npade
        self._theta = theta
        self._damped_wavenumber = damped_wavenumber
        self._parameters = parameters

    def _weak_form_impl(self):
        from bempp.api.operators.boundary.sparse import identity
        from bempp.api.operators.boundary.sparse import laplace_beltrami
        from bempp.api import ZeroBoundaryOperator
        from bempp.api import InverseSparseDiscreteBoundaryOperator

        import numpy as np

        space = self._space

        # Get the operators
        mass = identity(space, space, space,
                        parameters=self._parameters).weak_form()
        stiff = laplace_beltrami(
            space, space, space, parameters=self._parameters).weak_form()

        # Compute damped wavenumber

        if self._damped_wavenumber is None:
            bbox = self._space.grid.bounding_box
            rad = np.linalg.norm(bbox[1, :] - bbox[0, :]) / 2
            self._damped_wavenumber = self._wave_number + 1j * 0.4 * \
                self._wave_number ** (1.0 / 3.0) * rad ** (-2.0 / 3.0)

        # Get the Pade coefficients

        c0, alpha, beta = _pade_coeffs(self._npade, self._theta)

        # Now put all operators together

        term = ZeroBoundaryOperator(space, space, space).weak_form()

        for i in range(self._npade):
            term -= (alpha[i] / (self._damped_wavenumber ** 2) * stiff *
                     InverseSparseDiscreteBoundaryOperator(
                         mass - beta[i] / (self._damped_wavenumber ** 2) * stiff))
        result = 1j * self._wave_number * (c0 * mass + term * mass)

        return result


def osrc_ntd(space, wave_number, npade=2, theta=_np.pi / 3, damped_wavenumber=None, parameters=None, label="osrc_ntd"):
    """Return the OSRC approximation to the NtD map.

    Parameters
    ----------
    space : bempp.api.space.Space
        Underlying function space. Needs to be a continuous space.
    wave_number : complex
        Wavenumber of the operator.
    npade : int
        Number of Pade terms to be used in the approximation
    theta : float64
        Angle of the branch cut of the square root operator.
    damped_wavenumber : complex
        Damped wavenumber of the OSRC operator.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.
    label : string
        Label for the operator.

    """
    return _OsrcNtd(space, wave_number, npade, theta, damped_wavenumber, label=label)


class _OsrcNtd(_BoundaryOperator):

    def __init__(self, space, wave_number, npade=2, theta=_np.pi / 3., damped_wavenumber=None,
                 parameters=None, label="osrc_ntd"):
        super(_OsrcNtd, self).__init__(space, space, space, label=label)

        self._space = space
        self._wave_number = wave_number
        self._npade = npade
        self._theta = theta
        self._damped_wavenumber = damped_wavenumber
        self._parameters = parameters

    def _weak_form_impl(self):
        from bempp.api.operators.boundary.sparse import identity
        from bempp.api.operators.boundary.sparse import laplace_beltrami
        from bempp.api import ZeroBoundaryOperator
        from bempp.api import InverseSparseDiscreteBoundaryOperator

        import numpy as np

        space = self._space

        # Get the operators
        mass = identity(space, space, space,
                        parameters=self._parameters).weak_form()
        stiff = laplace_beltrami(
            space, space, space, parameters=self._parameters).weak_form()

        # Compute damped wavenumber

        if self._damped_wavenumber is None:
            bbox = self._space.grid.bounding_box
            rad = np.linalg.norm(bbox[1, :] - bbox[0, :]) / 2
            self._damped_wavenumber = self._wave_number + 1j * 0.4 * \
                self._wave_number ** (1.0 / 3.0) * rad ** (-2.0 / 3.0)

        # Get the Pade coefficients

        c0, alpha, beta = _pade_coeffs(self._npade, self._theta)

        # Now put all operators together

        term = ZeroBoundaryOperator(space, space, space).weak_form()

        for i in range(self._npade):
            term -= (alpha[i] / (self._damped_wavenumber ** 2) * stiff *
                     InverseSparseDiscreteBoundaryOperator(
                         mass - beta[i] / (self._damped_wavenumber ** 2) * stiff))
        result = 1. / (1j * self._wave_number) * (
            mass * InverseSparseDiscreteBoundaryOperator(mass - 1. / (self._damped_wavenumber ** 2) * stiff) * (c0 * mass + term * mass))

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

    import numpy as np

    # First define the coefficients for the standard Pade approximation.
    Np = n

    # C0 = 1.0
    aj = np.zeros(Np)
    bj = np.zeros(Np)
    for jj in range(1, Np + 1):  # remember Python skips the last index
        # remember that Python starts with zero in an array
        aj[jj - 1] = 2.0 / (2.0 * Np + 1.0) * np.sin(jj *
                                                     np.pi / (2.0 * Np + 1.0)) ** 2
        bj[jj - 1] = np.cos(jj * np.pi / (2.0 * Np + 1.0)) ** 2
    # Now define the coefficients for the 'theta' branch cut.
    C0t = np.exp(1j * theta / 2.0) * (
        1.0 + np.sum((aj * (np.exp(-1j * theta) - 1.0)) / (1.0 + bj * (np.exp(-1j * theta) - 1.0))))
    ajt = np.exp(-1j * theta / 2.0) * aj / \
        ((1 + bj * (np.exp(-1j * theta) - 1)) ** 2)
    bjt = np.exp(-1j * theta) * bj / (1 + bj * (np.exp(-1j * theta) - 1))
    return C0t, ajt, bjt
