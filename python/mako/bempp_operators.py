__all__ = ['dtypes', 'compatible_dtypes', 'bops']
from space import dtypes, compatible_dtypes

# Describes boundary operators
bops = {}
bops["laplace3dSingleLayerBoundaryOperator"] = {
    'location': ('laplace_3d', 'single_layer'),
    'doc': "Single-layer-potential boundary operator of 3D Laplace equation"
}
bops["laplace3dDoubleLayerBoundaryOperator"] = {
    'location': ('laplace_3d', 'double_layer'),
    'doc': "Double-layer-potential boundary operator of 3D Laplace equation"
}
bops["identityOperator"] = {
    'location': ('identity', ),
    'doc': "General identity operator"
}


for key, description in bops.iteritems():
    if 'implementation' not in description:
        description['implementation'] = 'standard'

    if 'header' not in description:
        f = lambda x: x if x.islower() else '_' + x.lower()
        description['header'] = 'bempp/assembly/%s.hpp' \
            % (key[0].lower() + ''.join([f(u) for u in key[1:]]))

    if 'c_creator' not in description:
        description['c_creator'] = "c_creator_" + key

    if 'py_creator' not in description:
        description['py_creator'] = "_create_" + key
