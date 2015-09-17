import bempp
import pickle

class TestPickleGlobalOptions(object):

    def test_pickle_global_options(self):

        serialized = pickle.dumps(bempp.global_parameters)
        restored = pickle.loads(serialized)
        orig = bempp.global_parameters
        self.__compare_options(orig, restored)

    def __compare_options(self, orig, restored):
        assert orig.assembly.boundary_operator_assembly_type \
                           == restored.assembly.boundary_operator_assembly_type
        assert orig.assembly.potential_operator_assembly_type \
                           == restored.assembly.potential_operator_assembly_type
        assert orig.assembly.enable_singular_integral_caching \
                           == restored.assembly.enable_singular_integral_caching

        assert orig.quadrature.double_singular \
                           == restored.quadrature.double_singular

        assert orig.quadrature.near.max_rel_dist \
                           == restored.quadrature.near.max_rel_dist
        assert orig.quadrature.near.single_order \
                           == restored.quadrature.near.single_order
        assert orig.quadrature.near.double_order \
                           == restored.quadrature.near.double_order

        assert orig.quadrature.medium.max_rel_dist \
                           == restored.quadrature.medium.max_rel_dist
        assert orig.quadrature.medium.single_order \
                           == restored.quadrature.medium.single_order
        assert orig.quadrature.medium.double_order \
                           == restored.quadrature.medium.double_order

        assert orig.quadrature.far.single_order \
                           == restored.quadrature.far.single_order
        assert orig.quadrature.far.double_order \
                           == restored.quadrature.far.double_order

        assert orig.hmat.assembly_mode \
                           == restored.hmat.assembly_mode
        assert orig.hmat.min_block_size \
                           == restored.hmat.min_block_size
        assert orig.hmat.max_block_size \
                           == restored.hmat.max_block_size
        assert orig.hmat.eta \
                           == restored.hmat.eta
        assert orig.hmat.eps \
                           == restored.hmat.eps
        assert orig.hmat.max_rank \
                           == restored.hmat.max_rank
        assert orig.hmat.compression_algorithm \
                           == restored.hmat.compression_algorithm
        assert orig.hmat.coarsening \
                           == restored.hmat.coarsening
        assert orig.hmat.coarseningAccuracy \
                           == restored.hmat.coarseningAccuracy
        assert orig.hmat.mat_vec_parallel_levels \
                           == restored.hmat.mat_vec_parallel_levels
