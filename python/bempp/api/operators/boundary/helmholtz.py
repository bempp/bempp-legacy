# pylint: disable-msg=too-many-arguments
from bempp.api.assembly.boundary_operator import BoundaryOperator as _BoundaryOperator
import numpy as _np

"""Definition of the Helmholtz boundary operators."""


def single_layer(domain, range_, dual_to_range,
                 wave_number,
                 label="SLP", symmetry='no_symmetry',
                 parameters=None):
    """Return the single-layer boundary operator."""

    from .modified_helmholtz import single_layer as sl

    return sl(domain, range_, dual_to_range,
              wave_number / (1j), label, symmetry,
              parameters)


def double_layer(domain, range_, dual_to_range,
                 wave_number,
                 label="DLP", symmetry='no_symmetry',
                 parameters=None):
    """Return the double-layer boundary operator."""

    from .modified_helmholtz import double_layer as dl

    return dl(domain, range_, dual_to_range,
              wave_number / (1j), label, symmetry,
              parameters)


def adjoint_double_layer(domain, range_, dual_to_range,
                         wave_number,
                         label="ADJ_DLP", symmetry='no_symmetry',
                         parameters=None):
    """Return the adjoint double-layer boundary operator."""

    from .modified_helmholtz import adjoint_double_layer as adl

    return adl(domain, range_, dual_to_range,
               wave_number / (1j), label, symmetry,
               parameters)


def hypersingular(domain, range_, dual_to_range,
                  wave_number,
                  label="HYP", symmetry='no_symmetry',
                  parameters=None,
                  use_slp=False):
    """Return the hypersingular boundary operator."""

    from .modified_helmholtz import hypersingular as hyp

    return hyp(domain, range_, dual_to_range,
               wave_number / (1j), label, symmetry,
               use_slp=use_slp, parameters=parameters)


def interior_calderon_projector(grid, wave_number, parameters=None):
    import bempp.api
    from bempp.api.assembly import BlockedOperator
    from bempp.api.space import project_operator

    if parameters is None:
        parameters = bempp.api.global_parameters

    blocked_operator = BlockedOperator(2, 2)

    const_space = bempp.api.function_space(grid, "DUAL", 0)
    lin_space = bempp.api.function_space(grid, "P", 1)
    lin_space_bary = bempp.api.function_space(grid, "B-P", 1)
    lin_space_disc_bary = bempp.api.function_space(grid, "B-DP", 1)
    lin_space_disc = bempp.api.function_space(grid.barycentric_grid(), "DP", 1)

    slp = bempp.api.operators.boundary.helmholtz.single_layer(lin_space_disc, lin_space_disc, lin_space_disc,
                                                            wave_number, parameters=parameters)
    ident1 = bempp.api.operators.boundary.sparse.identity(lin_space_bary, lin_space_bary, const_space)
    ident2 = bempp.api.operators.boundary.sparse.identity(const_space, const_space, lin_space_bary)
    dlp_disc = bempp.api.operators.boundary.helmholtz.double_layer(lin_space_disc, lin_space_disc, lin_space_disc,
                                                                 wave_number, parameters=parameters)
    dlp = project_operator(dlp_disc, domain=lin_space_bary, range_=lin_space_bary,
                           dual_to_range=const_space)

    adlp = project_operator(dlp_disc.transpose(const_space), domain=const_space, range_=const_space,
                            dual_to_range=lin_space_bary)

    blocked_operator[0, 1] = project_operator(slp, domain=const_space, range_=lin_space_bary,
                                              dual_to_range=const_space)

    blocked_operator[1, 0] = (bempp.api.operators.boundary.helmholtz.hypersingular( \
        lin_space_bary, const_space, lin_space_bary, wave_number,
        use_slp=project_operator(slp, domain=lin_space_disc_bary,
                                 dual_to_range=lin_space_disc_bary),
        parameters=parameters))
    blocked_operator[0, 0] = .5 * ident1 - dlp
    blocked_operator[1, 1] = .5 * ident2 + adlp

    return blocked_operator


def exterior_calderon_projector(grid, wave_number, parameters=None):
    import bempp.api
    from bempp.api.assembly import BlockedOperator
    from bempp.api.space import project_operator

    if parameters is None:
        parameters = bempp.api.global_parameters

    blocked_operator = BlockedOperator(2, 2)

    const_space = bempp.api.function_space(grid, "DUAL", 0)
    lin_space = bempp.api.function_space(grid, "P", 1)
    lin_space_bary = bempp.api.function_space(grid, "B-P", 1)
    lin_space_disc_bary = bempp.api.function_space(grid, "B-DP", 1)
    lin_space_disc = bempp.api.function_space(grid.barycentric_grid(), "DP", 1)

    slp = bempp.api.operators.boundary.helmholtz.single_layer(lin_space_disc, lin_space_disc, lin_space_disc,
                                                            wave_number, parameters=parameters)
    ident1 = bempp.api.operators.boundary.sparse.identity(lin_space_bary, lin_space_bary, const_space)
    ident2 = bempp.api.operators.boundary.sparse.identity(const_space, const_space, lin_space_bary)
    dlp_disc = bempp.api.operators.boundary.helmholtz.double_layer(lin_space_disc, lin_space_disc, lin_space_disc,
                                                                 wave_number, parameters=parameters)
    dlp = project_operator(dlp_disc, domain=lin_space_bary, range_=lin_space_bary,
                           dual_to_range=const_space)

    adlp = project_operator(dlp_disc.transpose(const_space), domain=const_space, range_=const_space,
                            dual_to_range=lin_space_bary)

    blocked_operator[0, 1] = -1. * project_operator(slp, domain=const_space, range_=lin_space_bary,
                                                    dual_to_range=const_space)

    blocked_operator[1, 0] = -1. * (bempp.api.operators.boundary.helmholtz.hypersingular( \
        lin_space_bary, const_space, lin_space_bary, wave_number
        use_slp=project_operator(slp, domain=lin_space_disc_bary,
                                 dual_to_range=lin_space_disc_bary),
        parameters=parameters))
    blocked_operator[0, 0] = .5 * ident1 + dlp
    blocked_operator[1, 1] = .5 * ident2 - adlp

    return blocked_operator



def osrc_dtn(space, wave_number, npade=2, theta=_np.pi / 3,
             parameters=None, label="osrc_dtn"):
    """Return the OSRC approximation to the DtN map."""
    return _OsrcDtn(space, wave_number, npade, theta, label=label)


class _OsrcDtn(_BoundaryOperator):
    def __init__(self, space, wave_number, npade=2, theta=_np.pi / 3.,
                 parameters=None, label="osrc_dtn"):
        super(_OsrcDtn, self).__init__(space, space, space, label=label)

        self._space = space
        self._wave_number = wave_number
        self._npade = npade
        self._theta = theta
        self._parameters = parameters

    def _weak_form_impl(self):
        from bempp.api.operators.boundary.sparse import identity
        from bempp.api.operators.boundary.sparse import laplace_beltrami
        from bempp.api import ZeroBoundaryOperator
        from bempp.api import InverseSparseDiscreteBoundaryOperator

        import numpy as np

        space = self._space

        # Get the operators
        mass = identity(space, space, space, parameters=self._parameters).weak_form()
        stiff = laplace_beltrami(space, space, space, parameters=self._parameters).weak_form()

        # Compute damped wavenumber

        bbox = self._space.grid.bounding_box
        rad = np.linalg.norm(bbox[1, :] - bbox[0, :]) / 2
        dk = self._wave_number + 1j * 0.4 * self._wave_number ** (1.0 / 3.0) * rad ** (-2.0 / 3.0)

        # Get the Pade coefficients

        c0, alpha, beta = _pade_coeffs(self._npade, self._theta)

        # Now put all operators together

        term = ZeroBoundaryOperator(space, space, space).weak_form()

        for i in range(self._npade):
            term -= (alpha[i] / (dk ** 2) * stiff *
                     InverseSparseDiscreteBoundaryOperator(
                         mass - beta[i] / (dk ** 2) * stiff))
        result = 1j * self._wave_number * (c0 * mass + term * mass)

        return result


def osrc_ntd(space, wave_number, npade=2, theta=_np.pi / 3, parameters=None, label="osrc_ntd"):
    """Return the OSRC approximation to the NtD map."""
    return _OsrcNtd(space, wave_number, npade, theta, label=label)


class _OsrcNtd(_BoundaryOperator):
    def __init__(self, space, wave_number, npade=2, theta=_np.pi / 3.,
                 parameters=None, label="osrc_ntd"):
        super(_OsrcNtd, self).__init__(space, space, space, label=label)

        self._space = space
        self._wave_number = wave_number
        self._npade = npade
        self._theta = theta
        self._parameters = parameters

    def _weak_form_impl(self):
        from bempp.api.operators.boundary.sparse import identity
        from bempp.api.operators.boundary.sparse import laplace_beltrami
        from bempp.api import ZeroBoundaryOperator
        from bempp.api import InverseSparseDiscreteBoundaryOperator

        import numpy as np

        space = self._space

        # Get the operators
        mass = identity(space, space, space, parameters=self._parameters).weak_form()
        stiff = laplace_beltrami(space, space, space, parameters=self._parameters).weak_form()

        # Compute damped wavenumber

        bbox = self._space.grid.bounding_box
        rad = np.linalg.norm(bbox[1, :] - bbox[0, :]) / 2
        dk = self._wave_number + 1j * 0.4 * self._wave_number ** (1.0 / 3.0) * rad ** (-2.0 / 3.0)

        # Get the Pade coefficients

        c0, alpha, beta = _pade_coeffs(self._npade, self._theta)

        # Now put all operators together

        term = ZeroBoundaryOperator(space, space, space).weak_form()

        for i in range(self._npade):
            term -= (alpha[i] / (dk ** 2) * stiff *
                     InverseSparseDiscreteBoundaryOperator(
                         mass - beta[i] / (dk ** 2) * stiff))
        result = 1. / (1j * self._wave_number) * (
            mass * InverseSparseDiscreteBoundaryOperator(mass - 1. / (dk ** 2) * stiff) * (c0 * mass + term * mass))

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
        aj[jj - 1] = 2.0 / (2.0 * Np + 1.0) * np.sin(jj * np.pi / (2.0 * Np + 1.0)) ** 2
        bj[jj - 1] = np.cos(jj * np.pi / (2.0 * Np + 1.0)) ** 2
    # Now define the coefficients for the 'theta' branch cut.
    C0t = np.exp(1j * theta / 2.0) * (
        1.0 + np.sum((aj * (np.exp(-1j * theta) - 1.0)) / (1.0 + bj * (np.exp(-1j * theta) - 1.0))))
    ajt = np.exp(-1j * theta / 2.0) * aj / ((1 + bj * (np.exp(-1j * theta) - 1)) ** 2)
    bjt = np.exp(-1j * theta) * bj / (1 + bj * (np.exp(-1j * theta) - 1))
    return C0t, ajt, bjt
