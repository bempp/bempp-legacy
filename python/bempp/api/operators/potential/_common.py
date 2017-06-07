"""Define common decorators for logging Potential operator assembly."""


def potential_logger(fun):
    """A decorator that inserts logging into the assembly of potentials."""
    import bempp.api
    import time

    def wrapper(*args, **kwargs):
        """Wrapper for potentials to emit log messages."""

        if 'parameters' in kwargs:
            parameters = kwargs['parameters']
        else:
            parameters = bempp.api.global_parameters

        mode = parameters.assembly.potential_operator_assembly_type
        bempp.api.log(
            ("POTENTIAL OPERATOR ASSEMBLY START. Number of points: %i." +
            "Space dimension: %i. Assembly type: %s ") % (
            len(args[1][0]), args[0].global_dof_count, mode))
        start = time.time()
        operator = fun(*args, **kwargs)
        end = time.time()
        bempp.api.log(
            "FINISHED POTENTIAL OPERATOR ASSEMBLY. Time: %.2E" % (
            end - start))
        return operator
    return wrapper
