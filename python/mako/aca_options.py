from collections import OrderedDict
properties = OrderedDict()
""" Attributes defined in AcaOptions """

properties['eps'] = 'double', 1e-4,\
    """ Estimate of the desired approximation accuracy """
properties['eta'] = 'double', 1e-4,\
    """ Cluster-pair admissibility parameter """
properties['minimumBlockSize'] = 'unsigned int', 16, \
    """ Minimum size of blocks approximated with ACA """
properties['maximumBlockSize'] = 'unsigned int', 'UINT_MAX-1', \
    """ Minimum size of blocks approximated with ACA """
properties['maximumRank'] = 'unsigned int', 'UINT_MAX-1', \
    """ Maximum rank of blocks stored in the low-rank format """
properties['recompress'] = 'bool', False, \
    """ Whether to recompress ACA matrix after construction """
properties['scaling'] = 'double', 1, \
    """ Estimate of the magnitude of the matrix elements """
properties['firstClusterIndex'] = 'int', -1, \
    """ Index of the first block cluster to approximate """

enums = {
    'AcaAssemblyMode': [
        "GLOBAL_ASSEMBLY", "LOCAL_ASSEMBLY", "HYBRID_ASSEMBLY"
    ],
    'ReactionToUnsupportedMode': ["IGNORE", "WARNING", "ERROR"]
}

enum_properties = {
    'AcaAssemblyMode': (
        'global', 'Discretization strategies during ACA assembly'
    ),
    'ReactionToUnsupportedMode': (
        'warning', 'Action when an unsupported assembly mode is detected.'
    )
}
