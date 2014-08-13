from bempp.fiber.quadrature_options cimport QuadratureOptions
def toggle_frozen(QuadratureOptions quadops):
    """ Toggles frozen status on and off """
    quadops.__is_frozen = not quadops.__is_frozen 
