import numpy as _np


class SoundHardPlaneWaveScattering(object):
    def __init__(self, space, wavenumber, coupling='osrc',
                 direction=_np.array([1., 0, 0])):
        self._space = space
        self._wavenumber = wavenumber
        self._coupling = coupling
        self._direction = direction
        self._id = None
        self._osrc = None
        self._double_layer = None
        self._hypersingular = None
        self._burton_miller_operator = None
        self._it_count = None

    @property
    def direction(self):
        return self._direction

    @direction.setter
    def direction(self, value):
        self._direction = value

    @property
    def wavenumber(self):
        return self._wavenumber

    @wavenumber.setter
    def wavenumber(self, value):

        self._wavenumber = value
        self._osrc = None
        self._double_layer = None
        self._hypersingular = None
        self._burton_miller_operator = None

    @property
    def iteration_count(self):
        return self._it_count

    @property
    def coupling(self):
        return self._coupling

    @coupling.setter
    def coupling(self, value):
        if value not in ['default', 'osrc']:
            raise ValueError("'coupling' must be one of: 'osrc', 'default'")
        self._coupling = value
        self._burton_miller_operator = None

    @property
    def space(self):
        return self._space

    @space.setter
    def space(self, value):

        self._space = space
        self._id = None
        self._double_layer = None
        self._hypersingular = None
        self._osrc = None
        self._burton_miller_operator = None

    def _compute_burton_miller(self):

        from bempp.operators.boundary.sparse import identity
        from bempp.operators.boundary.helmholtz import double_layer
        from bempp.operators.boundary.helmholtz import hypersingular
        from bempp.operators.boundary.helmholtz import osrc_ntd

        space = self._space
        wavenumber = self._wavenumber

        if self._id is None:
            self._id = identity(space, space, space)

        if self._double_layer is None:
            self._double_layer = double_layer(space, space,
                                              space, wavenumber)

        if self._hypersingular is None:
            self._hypersingular = hypersingular(space, space,
                                                space, wavenumber)

        if self._coupling == 'osrc' and self._osrc is None:
            self._osrc = osrc_ntd(space, wavenumber)

        if self._burton_miller_operator is None:
            if self._coupling == 'osrc':
                self._burton_miller_operator = (
                    .5 * self._id - self._double_layer - self._osrc * self._hypersingular)
            else:
                self._burton_miller_operator = (
                    .5 * self._id - self._double_layer - (1j / self._wavenumber) * self._hypersingular)

    def _compute_rhs(self):

        import bempp

        self._compute_burton_miller()

        def dirichlet_fun(x, n, domain_index, result):
            result[0] = _np.exp(1j * self._wavenumber * _np.dot(x, self._direction))

        def neumann_fun(x, n, domain_index, result):
            result[0] = _np.dot(n, 1j * self._wavenumber * self._direction * _np.exp(
                1j * self._wavenumber * _np.dot(x, self._direction)))

        g1 = bempp.GridFunction(self._space, fun=dirichlet_fun, complex_data=True)
        g2 = bempp.GridFunction(self._space, fun=neumann_fun, complex_data=True)

        if self._coupling == 'osrc':
            return -g1 + self._osrc * g2
        else:
            return -g1 + (1j / self._wavenumber) * g2

    def compute(self, tol=1E-5, maxiter=200):

        import scipy.sparse.linalg
        import bempp

        self._compute_burton_miller()
        self._it_count = 0

        def it_count(rk):
            self._it_count += 1

        x, _ = scipy.sparse.linalg.gmres(self._burton_miller_operator.strong_form(),
                                         self._compute_rhs().coefficients, tol=tol, maxiter=maxiter,
                                         callback=it_count)

        return bempp.GridFunction(self._space, coefficients=x)
