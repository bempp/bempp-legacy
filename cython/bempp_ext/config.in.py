""" Compiled time configurations  """
__all__ = ['paths', 'version']

version = "@Bempp_VERSION@"
""" BEM++ version string """


class paths:
    """ Holds paths to BEM++ installation """
    cmake_module = "@BEMPP_CMAKE_PATH@"
    """ Path with cmake module

        CMake projects should add this path to CMAKE_MODULE_PATH for CMake to
        find BEM++.
    """
    meshes = "@BEMPP_MESHES@"
    """ Path to bem++ example meshes

        Used in unit-tests and examples.
    """
