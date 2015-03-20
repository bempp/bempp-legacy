import dolfin as _dolfin

class FenicsOperator(object):

    def __init__(self,fenics_weak_form):

        self._fenics_weak_form = fenics_weak_form

    def weak_form(self):

        return _dolfin.assemble(self._fenics_weak_form).sparray()



