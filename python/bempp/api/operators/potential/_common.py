"""Define common decorators for logging Potential operator assembly."""

def potential_logger(fun):
    """A decorator that inserts logging into the assembly of potentials."""
    import bempp.api
    import time


    def wrapper(*args, **kwargs):

        if 'parameters' in kwargs:
            parameters = kwargs['parameters']
        else:
            parameters = bempp.api.global_parameters

        mode = parameters.assembly.potential_operator_assembly_type
        bempp.api.LOGGER.info("POTENTIAL OPERATOR ASSEMBLY START. Number of points: {0}. Space dimension: {1}. Assembly type: {2} ".format(
            len(args[1][0]), args[0].global_dof_count, mode))
        start = time.time()
        op = fun(*args, **kwargs)
        end = time.time()
        bempp.api.LOGGER.info("FINISHED POTENTIAL OPERATOR ASSEMBLY. Time: {0:.2E} sec.".format(end-start))
        return op
    return wrapper

