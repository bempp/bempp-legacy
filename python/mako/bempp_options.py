""" Defines assembly options """


def pythonic(name):
    """ Tries and creates pythonic name from c names """
    name = str(name)
    to_lower = lambda x: x if x == x.lower() else "_" + x.lower() 
    return name[0].lower() + ''.join([to_lower(u) for u in name[1:]])


# Descibes wrapper implementation
options = {}
options['SingularIntegralCaching'] = {
    'type': 'cbool',
    'default': True,
    'doc': 'Whether to cacher during weak form assembly',
}
options['SparseStorageOfLocalOperators'] = {
    'type': 'cbool',
    'default': True,
    'doc': 'Whether to store local operator in sparse format',
}
options['JointAssembly'] = {
    'type': 'cbool',
    'default': False,
    'doc': 'Integral operator superpositions are assembled jointly',
}
options['BlasInQuadrature'] = {
    'type': 'BlasQuadrature',
    'default': 'BLAS_QUADRATURE_AUTO',
    'doc': 'Use BLAS routines during quadrature',
    'doc_default': "'auto'",
    'getter': 'isBlasEnabledInQuadrature',
    'enums': {  # python name: cython name
        "auto": 'BLAS_QUADRATURE_AUTO',
        True: 'BLAS_QUADRATURE_YES',
        False: 'BLAS_QUADRATURE_NO'
    }
}
options['uniform_quadrature'] = {
    'type': 'cbool',
    'default': True,
    'doc': 'Homogeneous quadrature order over regular integrals in each block',
    'setter': 'makeQuadratureOrderUniformInEachCluster',
    'getter': 'isQuadratureOrderUniformInEachCluster'
}
options['Verbosity'] = {
    'type': 'VerbosityLevel',
    'default': 'VERBOSITY_MEDIUM',
    'doc': 'Verbosity during assembly',
    'doc_default': "'medium'",
    'setter': 'setVerbosityLevel',
    'getter': 'verbosityLevel',
    'enums': {  # python name: cython name
        "low": 'VERBOSITY_LOW',
        "medium": 'VERBOSITY_MEDIUM',
        "high": 'VERBOSITY_HIGH',
    }
}
options['mode'] = {
    'type': 'AcaAssemblyMode',
    'default': "GLOBAL_ASSEMBLY",
    'doc_default': "'global'",
    'enums': {
        "global": 'GLOBAL_ASSEMBLY',
        "local": 'LOCAL_ASSEMBLY',
        "hybrid": 'HYBRID_ASSEMBLY',
    },
    'pyname': 'aca_assembly_mode',
    'doc': 'Discretization strategies during ACA assembly'
}
options['reactionToUnsupportedMode'] = {
    'type': 'ReactionToUnsupportedMode',
    'default': "WARNING",
    'doc_default': "'warning'",
    'enums': {
        "warning": 'WARNING',
        "ignore": 'IGNORE',
        "error": 'ERROR',
    },
    'doc': 'Action when an unsupported assembly mode is detected'
}

options['eps'] = {
    'type': 'double',
    'default': 1e-4,
    'doc': 'Estimate of the desired approximation accuracy',
}
options['eta'] = {
    'type': 'double',
    'default': 1e-4,
    'doc': 'Cluster-pair admissibility parameter',
}
options['minimumBlockSize'] = {
    'type': 'unsigned int',
    'default': 16,
    'doc': 'Minimum size of blocks approximated with ACA',
}
options['maximumBlockSize'] = {
    'type': 'unsigned int',
    'default': 'UINT_MAX-1',
    'doc': 'Minimum size of blocks approximated with ACA',
}
options['maximumRank'] = {
    'type': 'unsigned int',
    'default': 'UINT_MAX-1',
    'doc': 'Maximum rank of blocks stored in the low-rank format',
}
options['recompress'] = {
    'type': 'bool',
    'default': False,
    'doc': 'Whether to recompress ACA matrix after construction',
}
options['scaling'] = {
    'type': 'double',
    'default': 1,
    'doc': 'Estimate of the magnitude of the matrix elements',
}
options['firstClusterIndex'] = {
    'type': 'int',
    'default': -1,
    'doc': 'Index of the first block cluster to approximate',
}
options['assembly_mode'] = {
    'doc_type': "'dense'|'aca'",
    'doc': "Switches between dense and aca assembly modes",
    'default': "'dense'",
    'implementation': 'manual',
    'c origin': 'manual'
}
options['max_threads'] = {
    'doc_type': "int|'auto'",
    'doc': 'Max number of threads during assembly',
    'default': "'auto'",
    'implementation': 'manual',
    'c origin': 'manual'
}
options['accuracy'] = {
    'doc_type': ":py:class:`AccuracyOptions`",
    'doc': 'Quadrature accuracy options',
    'default': "AccuracyOptions()",
    'implementation': 'manual',
    'c origin': 'manual'
}

# Add missing default values
assembly_options = [
    'Verbosity', 'uniform_quadrature', 'BlasInQuadrature',
    'SingularIntegralCaching', 'SparseStorageOfLocalOperators',
    'JointAssembly'
]
aca_ops = [
    'eps', 'eta', 'scaling', 'minimumBlockSize', 'maximumBlockSize',
    'maximumRank', 'recompress', 'scaling', 'firstClusterIndex'
]
for name in aca_ops:
    options[name]['implementation'] = 'getsetters'

for option, desc in options.iteritems():
    if 'enums' in desc and 'doc_type' not in desc:
        desc['doc_type'] = '|'.join([repr(u) for u in desc['enums'].keys()])
    if 'doc_type' not in desc:
        doc_type = 'bool' if desc['type'] == 'cbool' else desc['type']
        desc['doc_type'] = doc_type
    if 'doc_default' not in desc:
        desc['doc_default'] = desc['default']

    if 'c origin' not in desc:
        is_assembly = option in assembly_options
        desc['c origin'] = 'AssemblyOptions' if is_assembly else 'AcaOptions'

    if 'setter' not in desc:
        if desc['c origin'] == 'AcaOptions':
            desc['setter'] = option
        else:
            desc['setter'] = 'enable%s' % option

    if 'getter' not in desc:
        if desc['c origin'] == 'AcaOptions':
            desc['getter'] = option
        else:
            desc['getter'] = 'is%sEnabled' % option

    if 'enums' in desc and 'implementation' not in desc:
        desc['implementation'] = 'enums'
    if 'implementation' not in desc:
        desc['implementation'] = 'getsetters'
    if 'pyname' not in desc:
        desc['pyname'] = pythonic(option)
